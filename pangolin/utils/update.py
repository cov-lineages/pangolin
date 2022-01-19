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


def update(version_dictionary, data_dir=None):
    """
    Using the github releases API check for the latest current release
    of the set of depdencies provided e.g., pangolin, scorpio, pangolearn and
    constellations for complete --update and just pangolearn and constellations
    for --update_data.

    Dictionary keys must be one of pangolin, scorpio, pangolearn, or constellations

    Compare these to the currently running versions and if newer releases
    exist update to them accordingly (or do nothing if current).
    Afterwards, exit program safely with a 0 exit code.

    version_dictionary: dictionary keyed with dependency names and version for
                        that dependency
                        e.g.
    {pangolin: string containing the __version__ data for the currently
                      running pangolin module
    pangolearn: string containing the __version__ data for the imported
                       pangoLEARN data module
    scorpio: string containing the __version__ data for the imported
                       scorpio module
    constellations: string containing the __version__ data for the imported
                       constellations data module
    pango-designation: string containing the __version__ data for the imported
                       pango_designation data module}

    """
    package_names = {'pangolearn': 'pangoLEARN',
                     'pango-designation': 'pango_designation'
                    }

    # flag if any element is update if everything is the latest release
    # we want to just continue running
    for dependency, version in version_dictionary.items():

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
        latest_release_tidied = latest_release.strip('data release').lstrip('v').strip()
        latest_release_tidied = LooseVersion(latest_release_tidied)

        # clean up version numbers to remove leading v's and "data release"
        version = version.strip('data release').lstrip('v').strip()
        if dependency not in ['pangolin', 'scorpio', 'pango-designation',
                              'pangolearn', 'constellations']:
            raise ValueError("Dependency name for auto-update must be one "
                             "of: 'pangolin', 'pangolearn', scorpio', "
                             "'constellations', 'pango-designation'")

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
                    shutil.move(os.path.join(tempdir, extracted_dir, dependency_package), destination_directory)
            else:
                subprocess.run([sys.executable, '-m', 'pip', 'install', '--upgrade',
                                f"git+https://github.com/cov-lineages/{dependency}.git@{latest_release}"],
                                check=True,
                                stdout=subprocess.DEVNULL,
                                stderr=subprocess.DEVNULL)
            print(f"{dependency} updated to {latest_release}", file=sys.stderr)
        elif version > latest_release_tidied:
            print(f"{dependency} ({version}) is newer than latest stable "
                  f"release ({latest_release}), not updating.", file=sys.stderr)
        else:
            print(f"{dependency} already latest release ({latest_release})",
                    file=sys.stderr)
    sys.exit(0)

