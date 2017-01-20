#!/bin/bash

## BGC Tag for build
BGC_TAG=9847464

## Subdirectory in BGC repo to use
BGC_SUBDIR=.

## Available protocols for acquiring BGC source code
BGC_GIT_HTTP_ADDRESS=https://github.com/ACME-Climate/Ocean-BGC.git
BGC_GIT_SSH_ADDRESS=git@github.com:ACME-Climate/Ocean-BGC.git
BGC_SVN_ADDRESS=https://github.com/ACME-Climate/Ocean-BGC-src/tags
BGC_WEB_ADDRESS=https://github.com/ACME-Climate/Ocean-BGC-src/archive

GIT=`which git`
SVN=`which svn`
PROTOCOL=""

# BGC exists. Check to see if it is the correct version.
# Otherwise, flush the directory to ensure it's updated.
if [ -d BGC ]; then

	if [ -d .BGC_all/.git ]; then
		cd .BGC_all
		CURR_TAG=$(git rev-parse --short HEAD)
		cd ../
		if [ "${CURR_TAG}" == "${BGC_TAG}" ]; then
			echo "BGC version is current. Skip update"
		else
			unlink BGC
			rm -rf .BGC_all
		fi
	else
		unlink BGC
		rm -rf .BGC_all
	fi
fi

# BGC Doesn't exist, need to acquire souce code
# If might have been flushed from the above if, in the case where it was svn or wget that acquired the source.
if [ ! -d BGC ]; then
	if [ -d .BGC_all ]; then
		rm -rf .BGC_all
	fi

	if [ "${GIT}" != "" ]; then
		echo " ** Using git to acquire BGC source. ** "
		PROTOCOL="git ssh"
		git clone ${BGC_GIT_SSH_ADDRESS} .BGC_all &> /dev/null
		if [ -d .BGC_all ]; then
			cd .BGC_all
			git checkout ${BGC_TAG} &> /dev/null
			cd ../
			ln -sf .BGC_all/${BGC_SUBDIR} BGC
		else
			git clone ${BGC_GIT_HTTP_ADDRESS} .BGC_all &> /dev/null
			PROTOCOL="git http"
			if [ -d .BGC_all ]; then
				cd .BGC_all
				git checkout ${BGC_TAG} &> /dev/null
				cd ../
				ln -sf .BGC_all/${BGC_SUBDIR} BGC
			fi
		fi
	elif [ "${SVN}" != "" ]; then
		echo " ** Using svn to acquire BGC source. ** "
		PROTOCOL="svn"
		svn co ${BGC_SVN_ADDRESS}/${BGC_TAG} .BGC_all &> /dev/null
		ln -sf .BGC_all/${BGC_SUBDIR} BGC
	else
		echo " ** Using wget to acquire BGC source. ** "
		PROTOCOL="svn"
		BGC_ZIP_DIR=`echo ${BGC_TAG} | sed 's/v//g'`
		BGC_ZIP_DIR="BGC-src-${BGC_ZIP_DIR}"
		if [ ! -e .${BGC_TAG}.zip ]; then
			wget ${BGC_WEB_ADDRESS}/${BGC_TAG}.zip &> /dev/null
		fi
		unzip ${BGC_TAG}.zip &> /dev/null
		mv ${BGC_TAG}.zip .${BGC_TAG}.zip
		mv ${BGC_ZIP_DIR} .BGC_all
		ln -sf .BGC_all/${BGC_SUBDIR} BGC
	fi
fi

if [ ! -d BGC ]; then
	echo " ****************************************************** "
	echo " ERROR: Build failed to acquire BGC source."
	echo ""
	echo " Please ensure your proxy information is setup properly for"
	echo " the protocol you use to acquire BGC."
	echo ""
	echo " The automated script attempted to use: ${PROTOCOL}"
	echo ""
	if [ "${PROTOCOL}" == "git http" ]; then
		echo " This protocol requires setting up the http.proxy git config option."
	elif [ "${PROTOCOL}" == "git ssh" ]; then
		echo " This protocol requires having ssh-keys setup, and ssh access to git@github.com."
		echo " Please use 'ssh -vT git@github.com' to debug issues with ssh keys."
	elif [ "${PROTOCOL}" == "svn" ]; then
		echo " This protocol requires having svn proxys setup properly in ~/.subversion/servers."
	elif [ "${PROTOCOL}" == "wget" ]; then
		echo " This protocol requires having the http_proxy and https_proxy environment variables"
		echo " setup properly for your shell."
	fi
	echo ""
	echo " ****************************************************** "
fi
