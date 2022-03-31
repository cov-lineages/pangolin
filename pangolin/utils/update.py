#!/usr/bin/env python

import os
import sys
import json
import shutil
import tarfile
import subprocess
from urllib import request
from distutils.version import LooseVersion
from tempfile import TemporaryDirectory, tempdir

from pangolin.utils.log_colours import green,cyan,red

version_dict_keys = ['pangolin', 'scorpio', 'pangolin-data', 'constellations', 'pangolin-assignment']


def get_latest_release(dependency):
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
    latest_release_tarball = latest_release[0]['tarball_url']
    # extract and clean up latest release version
    latest_release = latest_release[0]['tag_name']
    return latest_release, latest_release_tarball


def pip_install_dep(dependency, release):
    """
    Use pip install to install a cov-lineages repository with the specificed release
    """
    url = f"git+https://github.com/cov-lineages/{dependency}.git@{release}"
    subprocess.run([sys.executable, '-m', 'pip', 'install', '--upgrade', url],
                   check=True,
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)


def install_pangolin_assignment():
    """
    If the pangolin-assignment repo has not been installed already then install the latest release.
    """
    try:
        import pangolin_assignment
        print(f"pangolin-assignment already installed with version {pangolin_assignment.__version__}; use --update or --update-data if you wish to update it.", file=sys.stderr)

    except:
        latest_release, tarball = get_latest_release('pangolin-assignment')
        pip_install_dep('pangolin-assignment', latest_release)
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
    for --update_data.  If pangolin-assignment has been added to the installation
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
                pip_install_dep(dependency, latest_release)
            print(f"{dependency} updated to {latest_release}", file=sys.stderr)
        elif version > latest_release_tidied:
            print(f"{dependency} ({version}) is newer than latest stable "
                  f"release ({latest_release}), not updating.", file=sys.stderr)
        else:
            print(f"{dependency} already latest release ({latest_release})",
                    file=sys.stderr)
    sys.exit(0)

