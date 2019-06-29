#!/usr/bin/env python

# stdlib
import os, sys, re
import glob
import json
import pprint
import argparse
import hashlib
import time
# since py3 split up urllib
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

# third party
import jinja2

"""
Renders Jinja2 templates using variables from dependency files

The behavior is not (yet) recursive.
"""

input_directory = "viral-ngs-template"
output_directory = "viral-ngs"
source_url = ""

dir_path = os.path.dirname(os.path.realpath(__file__))

class VersionString(object):
    """
        Class to validate and parse PEP440 version strings (also used by conda)
        Shortened and derived from: https://github.com/pypa/packaging/blob/16.7/packaging/version.py
    """

    VERSION_PATTERN = r"""
    (?P<prefix>v?)
    (?:
        (?:(?P<epoch>[0-9]+)!)?                           # epoch
        (?P<release>[0-9]+(?:\.[0-9]+)*)                  # release segment
        (?P<pre>                                          # pre-release
            [-_\.]?
            (?P<pre_l>(a|b|c|rc|alpha|beta|pre|preview))
            [-_\.]?
            (?P<pre_n>[0-9]+)?
        )?
        (?P<post>                                         # post release
            (?:-(?P<post_n1>[0-9]+))
            |
            (?:
                [-_\.]?
                (?P<post_l>post|rev|r)
                [-_\.]?
                (?P<post_n2>[0-9]+)?
            )
        )?
        (?P<dev>                                          # dev release
            [-_\.]?
            (?P<dev_l>dev)
            [-_\.]?
            (?P<dev_n>[0-9]+)?
        )?
    )
    (?:\+(?P<local>[a-z0-9]+(?:[-_\.][a-z0-9]+)*))?       # local version
    """
    version_re = re.compile(
        r"^\s*" + VERSION_PATTERN + r"\s*$",
        re.VERBOSE | re.IGNORECASE,)

    def __init__(self, v):
        self.v = v

    def __str__(self):
        parts = []

        try:
            # 'v' prefix
            if self.version_re.match(self.v).group("prefix") is not None:
                parts.append("{0}".format(self.version_re.match(self.v).group("prefix")))

            # Epoch
            if ( int(self.version_re.match(self.v).group("epoch")) if self.version_re.match(self.v).group("epoch") else 0) != 0:
                parts.append("{0}!".format(self.version_re.match(self.v).group("epoch")))

            # Release segment
            parts.append(".".join(str(x) for x in self.version_re.match(self.v).group("release").split(".")))

            # Pre-release
            if self.version_re.match(self.v).group("pre") is not None:
                parts.append("".join(str(x) for x in self.version_re.match(self.v).group("pre")))

            # Post-release
            if self.version_re.match(self.v).group("post") is not None:
                parts.append("{0}".format(self.version_re.match(self.v).group("post")))

            # Development release
            if self.version_re.match(self.v).group("dev") is not None:
                parts.append("{0}".format(self.version_re.match(self.v).group("dev")))

            # Local version segment
            if self.version_re.match(self.v).group("local") is not None:
                parts.append(
                    "+{0}".format("".join(str(x) for x in self.version_re.match(self.v).group("local")))
                )
        except:
            raise argparse.ArgumentTypeError("String '%s' does not match required PEP440 format" % (self.v,))

        return "".join(parts)


def reformat_package_line(line):
    """
    This function is meant to take a package spec in conda or pip format
    and return one in conda recipe format: https://conda.io/docs/spec.html
    """
    # regex to match comment-only line
    comment_re = re.compile(r"^(?:\s*\#.*)$")

    # regex to match package spec line, with support for comments and selectors.
    # This will also capture hash-indicated selectors and comments (ex. "# [osx]")
    # which may, or may not, be useful in their original context.
    package_re = re.compile(r"^(?P<package>[a-zA-Z0-9\-\_]+)(?:\s*)(?:(?P<comparator>[\>\<=]?=?)(?:\s*)(?P<version>[^\s\#=]+)(?:=(?P<build>[0-9]*))?(?:\s*))?(?P<selector>\s*\#\s*\[.*\])?(?P<comment>\s*\#.*)?$")

    # when we need to specify a different comparator for the recipe
    comparator_replacements = {
        "=": "==",
    }

    # the line shold not have a newline
    line = line.replace("\n","").replace("\r","")

    # if the line is a comment, simpy return it
    if len(line)==0 or comment_re.match(line):
        return line
    # otherwise, build a package spec string suitable for a conda recipe
    else: 
        m = package_re.match(line)
        recipe_package_string = "- {package} {comparator}{version}{build}{selector}{comment}".format(
            package    = m.group("package").lower(), # conda packages must have lowercase names
            comparator = "" if not m.group("comparator") else comparator_replacements.get(m.group("comparator"), m.group("comparator")),
            version    = "" if not m.group("version") else m.group("version"),
            build      = "" if not m.group("build") else " "+m.group("build")+"*", # Todo: verify build separator character for recip format
            selector   = "" if not m.group("selector") else " "+m.group("selector"),
            comment    = "" if not m.group("comment") else " "+m.group("comment")
        )
        return recipe_package_string

def url_md5(url):
    hash_md5 = hashlib.md5()
    CHUNK_SIZE = 16 * 1024

    # try four times to download the file. If one fails, wait two seconds and try again.
    try_count = 1
    while True:
        try:
            print("Downloading source package for hash calculation...")
            print(url)
            response = urlopen(url)
            for chunk in iter(lambda: response.read(CHUNK_SIZE), b""):
                hash_md5.update(chunk)
            break
        except:
            print("Download {} failed, sleeping then retrying...".format(try_count))
            try_count +=1
            if try_count >3:
                raise
            time.sleep(2)
            continue

    return hash_md5.hexdigest()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Renger the conda recipe.')
    parser.add_argument('version',
                        type=VersionString,
                        help='the version number of the package')
    parser.add_argument('--package-name', dest="package_name",
                        type=str,
                        help='override the default name of the package described in the recipe')
    parser.add_argument('--build-reqs', nargs='*', dest='build_requirements',
                        type=argparse.FileType('r'),
                        help='build-time requirements file')
    parser.add_argument('--run-reqs', nargs='*', dest='run_requirements',
                        type=argparse.FileType('r'),
                        help='run-time requirements file')
    parser.add_argument('--py2-run-reqs', nargs='*', dest='py2_run_requirements',
                        type=argparse.FileType('r'),
                        help='python2-only run-time requirements file')
    parser.add_argument('--py3-run-reqs', nargs='*', dest='py3_run_requirements',
                        type=argparse.FileType('r'),
                        help='python3-only run-time requirements file')
    parser.add_argument('--linux-run-reqs', nargs='*', dest='linux_run_requirements',
                        type=argparse.FileType('r'),
                        help='linux-only run-time requirements file')
    parser.add_argument('--osx-run-reqs', nargs='*', dest='osx_run_requirements',
                        type=argparse.FileType('r'),
                        help='osx-only run-time requirements file')
    parser.add_argument('--test-reqs', nargs='*', dest='test_requirements',
                        type=argparse.FileType('r'),
                        help='test-time requirements file')
    parser.add_argument('--download-filename', dest='src_download_filename',
                        type=str,
                        help='An argument to override the usual filename to download; '
                        'useful for specifying a branch name.')

    try:
       args = parser.parse_args()
       if not any(vars(args).values()):
            parser.print_help()
            sys.exit(0)
    except:
        sys.exit(0)

    args_dict = vars(args)

    recipe_variables = {}

    # store two separate version strings, one to use for the conda package and one
    # that should match github tagged releases
    recipe_variables["PKG_VERSION"] = str(args_dict.pop("version"))
    if "src_download_filename" in args_dict and args_dict["src_download_filename"] is not None:
        recipe_variables["SRC_FILE_PREFIX"] = str(args_dict.pop("src_download_filename"))
    else:
        # otherwise use the "tag"
        recipe_variables["SRC_FILE_PREFIX"] = recipe_variables["PKG_VERSION"]

    if "package_name" in args_dict and args_dict["package_name"] is not None:
        recipe_variables["PKG_NAME"] = str(args_dict.pop("package_name"))

    # strip "v" prefix from versions that look like v1.14.0
    if recipe_variables["PKG_VERSION"].startswith("v"):
        recipe_variables["PKG_VERSION_CONDA"] = recipe_variables["PKG_VERSION"][1:]
    else:
        recipe_variables["PKG_VERSION_CONDA"] = recipe_variables["PKG_VERSION"]

    # after we pop the positional argument(s), the optional ones remaining are all files
    for var_name, req_files in args_dict.items():
        if req_files:
            for reqs_file in req_files:
                if reqs_file:
                    recipe_variables[var_name] = []
                    for line in reqs_file:
                        conda_style_package_line = reformat_package_line(line)
                        if len(conda_style_package_line):
                            recipe_variables[var_name].append(conda_style_package_line)
    pprint.pprint(recipe_variables)

    j_env = jinja2.Environment(loader=jinja2.FileSystemLoader(os.path.join(dir_path, input_directory)))

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    template_files = os.listdir(os.path.join(dir_path,input_directory))

    for template_file in template_files:
        print("Rendering "+ template_file)
        # jinja expects the filename to be just that, not a path
        # it should be relative to the FileSystemLoader() path set above
        template = j_env.get_template(template_file)
        output_from_parsed_template = template.render(recipe_variables)

        # save the rendered output
        with open(os.path.join(dir_path, output_directory, template_file), "w") as f:
            f.write(output_from_parsed_template)

        # populate md5 hashes for any source urls present
        if(template_file.endswith(".yaml")):
            # calculate and add md5 hashes to recipe
            with open(os.path.join(dir_path, output_directory, template_file), "r") as inf:
                with open(os.path.join(dir_path, output_directory, template_file+".checksumed"), "w") as outf:
                    for line in inf:
                        # if this is an md5 line, don't write it out
                        if line.strip().startswith("md5"):
                            continue
                        # if this is not an md5 line, write it verbatim
                        else:
                            outf.writelines([line])

                            # if this is a url line
                            if line.strip().startswith("url"):
                                # parse out the url
                                url_re = re.compile(r"^(?:(?P<leadingspace>\s*)url:\s*)(?P<url>[\S]*)(?P<extra>.*)$")
                                matches = url_re.match(line)
                                if matches:
                                    if matches.group("url"):
                                        # download file and calculate md5
                                        src_hash = url_md5(matches.group("url"))
                                        hash_line = "{leadingspace}md5: {src_hash}{extra}".format(
                                            leadingspace="" if not matches.group("leadingspace") else matches.group("leadingspace"),
                                            src_hash=src_hash,
                                            extra="" if not matches.group("extra") else matches.group("extra")
                                        )
                                        outf.writelines([hash_line+"\n"])

                                    else:
                                        raise Exception("The yaml file url line does not appear to contain a url")


            # move the file with checksums
            os.rename(os.path.join(dir_path, output_directory, template_file+".checksumed"), os.path.join(dir_path, output_directory, template_file))




