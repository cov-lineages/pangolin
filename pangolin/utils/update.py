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


def git_lfs_install():
    """
    'git-lfs install' must be run after installing git-lfs and before cloning a repo
    that uses Git LFS.
    """
    try:
        subprocess.run(['git-lfs', 'install'],
                   check=True,
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        stderr = e.stderr.decode('utf-8')
        sys.stderr.write(cyan(f"Error: {e}:\n{stderr}\n"))
        sys.exit(-1)

def pip_install_dep(dependency, release, datadir=None):
    """
    Use pip install to install a cov-lineages repository with the specificed release
    """
    url = f"git+https://github.com/cov-lineages/{dependency}.git@{release}"
    pip_command = [sys.executable, '-m', 'pip', 'install', '--upgrade']
    if datadir is not None:
        pip_command.extend(['--target', datadir])
    pip_command.append(url)
    subprocess.run(pip_command,
                   check=True,
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)


def install_pangolin_assignment(pangolin_assignment_version, datadir=None):
    """
    If the pangolin-assignment repo has not been installed already then install the latest release.
    """
    if pangolin_assignment_version is not None:
        print(f"pangolin-assignment already installed with version {pangolin_assignment_version}; use --update or --update-data if you wish to update it.", file=sys.stderr)
    else:
        git_lfs_install()
        latest_release, tarball = get_latest_release('pangolin-assignment')
        pip_install_dep('pangolin-assignment', latest_release, datadir)
        print(f"pangolin-assignment installed with latest release ({latest_release})")

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
            pip_install_dep(dependency, latest_release, data_dir)
            print(f"{dependency} updated to {latest_release}", file=sys.stderr)
        elif version > latest_release_tidied:
            print(f"{dependency} ({version}) is newer than latest stable "
                  f"release ({latest_release}), not updating.", file=sys.stderr)
        else:
            print(f"{dependency} already latest release ({latest_release})",
                    file=sys.stderr)
    sys.exit(0)

