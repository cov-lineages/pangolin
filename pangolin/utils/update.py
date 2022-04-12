#!/usr/bin/env python

import os
import sys
import json
import re
import shutil
import tarfile
import subprocess
from urllib import request
from distutils.version import LooseVersion
from tempfile import TemporaryDirectory, tempdir

from pangolin.utils.log_colours import green,cyan,red

version_dict_keys = ['pangolin', 'scorpio', 'pangolin-data', 'constellations', 'pangolin-assignment']

dependency_web_dir = { 'pangolin-assignment': 'https://hgdownload-test.gi.ucsc.edu/goldenPath/wuhCor1/pangolin-assignment' }


def get_latest_cov_lineages(dependency):
    """
    Using the github releases API check for the latest release of dependency and its tarball
    """
    try:
        latest_release = request.urlopen(
            f"https://api.github.com/repos/cov-lineages/{dependency}/releases")
    # to catch and give a useful error message when people try to run this
    # either update option on clusters without external connectivity
    # or have exceeded the github API limit temporarily
    # this may also catch genuine bugs when version and release tags diverge
    # so if this is thrown and there is definitely connectivity then
    # double check the version labels
    except Exception as e:
        sys.stderr.write(cyan("Unable to connect to reach github API "
                               "--update/--data_update requires internet "
                               "connectivity so may not work on certain "
                               "systems or if your IP has exceeded the "
                              f"5,000 request per hour limit\n{e}\n"))
        sys.exit(-1)

    latest_release = json.load(latest_release)
    try:
        # Find the latest stable release
        latest_release_dict = next(x for x in latest_release if not x['draft'] and not x['prerelease'])
    except:
        # All releases to date are prerelease or draft, just take the latest
        latest_release_dict = latest_release[0]
    latest_release_tarball = latest_release_dict['tarball_url']
    # extract and clean up latest release version
    latest_release = latest_release_dict['tag_name']
    return latest_release, latest_release_tarball


def get_latest_web_dir(dependency, web_dir):
    """
    Find the tarball url with the latest release from a web directory with versioned tarballs
    instead of github.  An HTTP GET of the web directory must return some text that contains
    names of files in that directory, some of which are {dependency}-{version}.tar.gz.
    """
    try:
        listing = request.urlopen(web_dir).read().decode('utf-8')
    except:
        sys.stderr.write(cyan(f"Unable to read {web_dir}"))
        sys.exit(-1)
    tarRe = re.compile(f"{dependency}-(.*?).tar.gz")
    matches = list(set(tarRe.findall(listing)))
    if not matches:
        sys.stderr.write(cyan(f"Can't find {dependency}-<version>.tar.gz files in listing of {web_dir}"))
        sys.exit(-1)
    versions = [LooseVersion(v) for v in matches]
    versions.sort()
    latest_release = str(versions[-1])
    latest_release_tarball = f"{web_dir}/{dependency}-{latest_release}.tar.gz"
    return latest_release, latest_release_tarball


def get_latest_release(dependency):
    """
    If dependency comes from a web directory then find latest release and tarball there, otherwise
    query github API for cov-lineages repo
    """
    if dependency in dependency_web_dir:
        return get_latest_web_dir(dependency, dependency_web_dir[dependency])
    else:
        return get_latest_cov_lineages(dependency)


def pip_install_url(url):
    """
    Use pip install to install a package from a url.
    """
    subprocess.run([sys.executable, '-m', 'pip', 'install', '--upgrade', url],
                   check=True,
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)


def pip_install_cov_lineages(dependency, release):
    """
    Use pip install to install a cov-lineages repository with the specified release
    """
    url = f"git+https://github.com/cov-lineages/{dependency}.git@{release}"
    pip_install_url(url)


def install_pangolin_assignment():
    """
    If the pangolin-assignment repo has not been installed already then install the latest release.
    """
    try:
        import pangolin_assignment
        print(f"pangolin-assignment already installed with version {pangolin_assignment.__version__}; use --update or --update-data if you wish to update it.", file=sys.stderr)

    except:
        latest_release, tarball = get_latest_release('pangolin-assignment')
        pip_install_url(tarball)
        print(f"pangolin-assignment installed with latest release ({latest_release})")


def add_pangolin_assignment_if_installed(version_dictionary):
    """
    If pangolin_assignment has been installed then add it to version_dictionary, else ignore.
    """
    try:
        import pangolin_assignment
        version_dictionary["pangolin-assignment"] = pangolin_assignment.__version__
    except:
        pass


def update(version_dictionary, data_dir=None):
    """
    Using the github releases API check for the latest current release
    of the set of dependencies provided e.g., pangolin, scorpio, pangolin-data and
    constellations for complete --update and just pangolearn and constellations
    for --update_data.  If pangolin-assignment has been added to version_dictionary
    then it will be included in both --update and --update-data.

    Dictionary keys must be one of pangolin, scorpio, pangolin-data, constellations
    or pangolin-assignment

    Compare these to the currently running versions and if newer releases
    exist update to them accordingly (or do nothing if current).
    Afterwards, exit program safely with a 0 exit code.

    version_dictionary: dictionary keyed with dependency names and version for
                        that dependency
                        e.g.
    {pangolin: string containing the __version__ data for the currently
                      running pangolin module
    pangolin-data: string containing the __version__ data for the imported
                       pangolin-data module
    scorpio: string containing the __version__ data for the imported
                       scorpio module
    constellations: string containing the __version__ data for the imported
                       constellations data module}

    """
    package_names = {'pangolin-data': 'pangolin_data',
                     'pangolin-assignment': 'pangolin_assignment'
                    }

    # flag if any element is update if everything is the latest release
    # we want to just continue running
    for dependency, version in version_dictionary.items():

        latest_release, latest_release_tarball = get_latest_release(dependency)
        latest_release_tidied = latest_release.strip('data release').lstrip('v').strip()
        latest_release_tidied = LooseVersion(latest_release_tidied)

        # clean up version numbers to remove leading v's and "data release"
        version = version.strip('data release').lstrip('v').strip()
        if dependency not in version_dict_keys:
            
            raise ValueError("Dependency name for auto-update must be one of: '" +
                             "', '".join(version_dict_keys) + "'")

        # convert to LooseVersion to have proper ordering of versions
        # this prevents someone using the latest commit/HEAD from being
        # downgraded to the last stable release
        version = LooseVersion(version)

        if version < latest_release_tidied:
            if data_dir is not None:
                # this path only gets followed when the user has --update_data and they
                # have also specified a --datadir
                with TemporaryDirectory() as tempdir:
                    dependency_package = package_names.get(dependency, dependency)
                    tarball_path = os.path.join(tempdir, 'tarball.tgz')
                    open(tarball_path, 'wb').write(request.urlopen(latest_release_tarball).read())
                    tf = tarfile.open(tarball_path)
                    extracted_dir = tf.next().name
                    tf.extractall(path=tempdir)
                    tf.close()
                    destination_directory = os.path.join(data_dir, dependency_package)
                    if os.path.isdir(destination_directory):
                        shutil.rmtree(destination_directory)
                    shutil.move(os.path.join(tempdir, extracted_dir, dependency_package), destination_directory)
            else:
                if dependency in dependency_web_dir:
                    pip_install_url(latest_release_tarball)
                else:
                    pip_install_cov_lineages(dependency, latest_release)
            print(f"{dependency} updated to {latest_release}", file=sys.stderr)
        elif version > latest_release_tidied:
            print(f"{dependency} ({version}) is newer than latest stable "
                  f"release ({latest_release}), not updating.", file=sys.stderr)
        else:
            print(f"{dependency} already latest release ({latest_release})",
                    file=sys.stderr)
    sys.exit(0)

