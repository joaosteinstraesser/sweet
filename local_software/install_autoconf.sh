#! /bin/bash

source ./install_helpers.sh ""

# Name of package
PKG_NAME="autoconf"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/autoconf"

# URL to source code to fetch it
PKG_URL_SRC="autoconf-2.69.tar.gz"

if [ "`uname -s`" != "Linux" ] && [ "`uname -s`" != "Darwin" ]; then
	echo "This script only supports Autoconf on Linux systems"
	exit 1
fi

config_setup

config_package $@

config_configure_make_default_install

config_success
