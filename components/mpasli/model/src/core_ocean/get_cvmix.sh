#!/bin/bash

## CVMix Tag for build
CVMIX_TAG=v0.64-beta

## Subdirectory in CVMix repo to use
CVMIX_SUBDIR=src/shared

## Available protocols for acquiring CVMix source code
CVMIX_GIT_HTTP_ADDRESS=https://github.com/CVMix/CVMix-src.git
CVMIX_GIT_SSH_ADDRESS=git@github.com:CVMix/CVMix-src.git
CVMIX_SVN_ADDRESS=https://github.com/CVMix/CVMix-src/tags
CVMIX_WEB_ADDRESS=https://github.com/CVMix/CVMix-src/archive

GIT=`which git`
SVN=`which svn`
PROTOCOL=""

# CVMix exists. Need to make sure it's updated if it is git.
# Otherwise, flush the directory to ensure it's updated.
if [ -d cvmix ]; then
	unlink cvmix

	if [ -d .cvmix_all/.git ]; then
		cd .cvmix_all
		git fetch origin &> /dev/null
		git checkout ${CVMIX_TAG} &> /dev/null
		cd ../
		ln -sf .cvmix_all/${CVMIX_SUBDIR} cvmix
	else
		rm -rf .cvmix_all
	fi
fi

# CVmix Doesn't exist, need to acquire souce code
# If might have been flushed from the above if, in the case where it was svn or wget that acquired the source.
if [ ! -d cvmix ]; then 
	if [ -d .cvmix_all ]; then
		rm -rf .cvmix_all
	fi

	if [ "${GIT}" != "" ]; then 
		echo " ** Using git to acquire cvmix source. ** "
		PROTOCOL="git https"
		git clone ${CVMIX_GIT_HTTP_ADDRESS} .cvmix_all &> /dev/null
		if [ -d .cvmix_all ]; then 
			cd .cvmix_all 
			git checkout ${CVMIX_TAG} &> /dev/null
			cd ../ 
			ln -sf .cvmix_all/${CVMIX_SUBDIR} cvmix 
		else 
			git clone ${CVMIX_GIT_SSH_ADDRESS} .cvmix_all &> /dev/null
			PROTOCOL="git ssh"
			if [ -d .cvmix_all ]; then 
				cd .cvmix_all 
				git checkout ${CVMIX_TAG} &> /dev/null
				cd ../ 
				ln -sf .cvmix_all/${CVMIX_SUBDIR} cvmix 
			fi 
		fi 
	elif [ "${SVN}" != "" ]; then 
		echo " ** Using svn to acquire cvmix source. ** "
		PROTOCOL="svn"
		svn co ${CVMIX_SVN_ADDRESS}/${CVMIX_TAG} .cvmix_all &> /dev/null
		ln -sf .cvmix_all/${CVMIX_SUBDIR} cvmix
	else 
		echo " ** Using wget to acquire cvmix source. ** "
		PROTOCOL="svn"
		CVMIX_ZIP_DIR=`echo ${CVMIX_TAG} | sed 's/v//g'`
		CVMIX_ZIP_DIR="CVMix-src-${CVMIX_ZIP_DIR}"
		if [ ! -e .${CVMIX_TAG}.zip ]; then
			wget ${CVMIX_WEB_ADDRESS}/${CVMIX_TAG}.zip &> /dev/null
		fi
		unzip ${CVMIX_TAG}.zip &> /dev/null
		mv ${CVMIX_TAG}.zip .${CVMIX_TAG}.zip
		mv ${CVMIX_ZIP_DIR} .cvmix_all 
		ln -sf .cvmix_all/${CVMIX_SUBDIR} cvmix 
	fi 
fi

if [ ! -d cvmix ]; then
	echo " ****************************************************** "
	echo " ERROR: Build failed to acquire CVMix source."
	echo ""
	echo " Please ensure your proxy information is setup properly for"
	echo " the protocol you use to acquire CVMix."
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
