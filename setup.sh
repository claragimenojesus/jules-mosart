#!/bin/bash
# Clara Gimeno Jésus, June 2022

## Bash script with instructions for setting up environments and directories before running JULES

# Sign up and the JULES source code via https://jules-lsm.github.io/access_req/JULES_access.html
# It can take a while (a week or two)

# Follow 'JULES from scratch' manual to install Cylc, Rose and FCM and cache your morsr password

## WB: seems like we need python2, so let's create a virtualenv first:

wget https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
bash Anaconda3-2022.05-Linux-x86_64.sh
conda install -c conda-forge pygraphviz # this solved some error messages

# make sure you copy the file environment.yml from the git repo
# note that we have two environments, one for jules and one for mosart!

conda env create -f environment-jules.yml

# Once donwload complete, run
conda activate jules

# create a directory for local software:

mkdir ~/local
cd ~/local

# Install Cylc    ## WB: can we use package cylc-flow?

git clone https://github.com/cylc/cylc.git

cd cylc
#git tag -l

git checkout tags/6.0.2
cd ..
export PATH=$HOME/.local/cylc/bin:$PATH
cylc --version
cd ~/local/cylc
cd ../..

export PATH=$HOME/.local/bin:$PATH
cylc check-software

## Install Rose
cd ~/.local
git clone git@github.com:metomi/rose.git

cd rose
#git tag -l
git checkout tags/2018.02.0     # Clara and Simon use more recent version (2019.01)
cd ..
export PATH=$HOME/.local/rose/bin:$PATH
rose --version
cd ..
mkdir ~/.metomi

# Create a text file ~/.metomi/rose.conf containing the following test and substituting your username and whatever your SITE is
echo "[rosie-id]
prefix-username.u=yourusername

[rose stem]
automatic-options=SITE=jasmin" > rose.conf.txt

cd ~/.local/rose/etc/

echo "[rosie-id]
prefix-default=u
prefixes-ws-default=u
prefix-location.u=https://code.metoffice.gov.uk/svn/roses-u
prefix-web.u=https://code.metoffice.gov.uk/trac/roses-u/intertrac/source:
prefix-ws.u=https://code.metoffice.gov.uk/rosie/u" >> rose.conf


.~/.bashrc

# Check the Rose installation and server links
cd ~/.local/
git clone git@github.com:metomi/metomi-vms.git
export PATH=$PATH:$HOME/.local/metomi-vms/usr/local/bin     #WB: may need to be installed. Note: at the end because it also comes with cylc
# note: you may need to enable gpg-preset-passphrase in .gnupg/gpg-agent.conf 
mosrs-cache-password
rosie hello                                                 

## Install FCM
## WB: the following code does not install fcm. Can be installed via apt install fcm
cd ~/.local
git clone https://github.com/metomi/fcm.git
cd fcm
git tag -l
git checkout tags/2017.10.0
cd ..
export PATH=$HOME/.local/fcm/bin:$PATH
fcm --version

ls ~/.subversion/servers
mkdir ~/.subversion
echo "[groups] 
metofficesharedrepos =code*.metoffice.gov.uk
[metofficesharedrepos]
username =yourusername
store-plaintext-passwords=no" >> servers.txt

# Check FCM installation
fcm --version

# Set up keywords
cd ~/.metomi/fcm/
echo "location{primary, type:svn}[jules.x] = https://code.metoffice.gov.uk/svn/jules/main
browser.loc-tmpl[jules.x] = https://code.metoffice.gov.uk/trac/{1}/intertrac/source:/{2}{3}
browser.comp-pat[jules.x] = (?msx-i:\A // [^/]+ /svn/ ([^/]+) /*(.*) \z)

location{primary, type:svn}[jules_doc.x] = https://code.metoffice.gov.uk/svn/jules/doc
browser.loc-tmpl[jules_doc.x] = https://code.metoffice.gov.uk/trac/{1}/intertrac/source:/{2}{3}
browser.comp-pat[jules_doc.x] = (?msx-i:\A // [^/]+ /svn/ ([^/]+) /*(.*) \z)" >> keyword.cfg 

## Caching your MOSRS password
#Caching your MOSRS password may be problematic if your MOSRS password has special characters in it (eg "%","!" or "&")
#Download two MOSRS utilities from [here]# (https://code.metoffice.gov.uk/trac/home/wiki/AuthenticationCaching/GpgAgent):
#firefox https://code.metoffice.gov.uk/trac/home/raw-attachment/wiki/AuthenticationCaching/GpgAgent/mosrs-setup-gpg-agent
#Put in your MOSRS password, save the file and close Firefox.
#firefox https://code.metoffice.gov.uk/trac/home/raw-attachment/wiki/AuthenticationCaching/GpgAgent/mosrs-cache-password
#Again, put in your MOSRS password, save the file and close Firefox.

mv ~/Downloads/mosrs-setup-gpg-agent ~/.local/mosrs-setup-gpg-agent

mv ~/Downloads/mosrs-cache-password ~/.local/mosrs-cache-password

export PATH=$HOME/.local:$PATH

chmod -R 755 ~/.local

#Now activate
. mosrs-setup-gpg-agent

#To test this has worked correctly
rosie hello
#svn info --username=[yourusername] https://code.metoffice.gov.uk/svn/test

#you will need to re-cache your MOSRS password every time you log in from now on.
export PATH=$HOME/.local:$PATH
[[ "$-" != *i* ]] && return
. mosrs-setup-gpg-agent


## Install JULES
mkdir ~/MODELS
cd ~/MODELS

# Download a version of JULES

fcm co fcm:jules.x_tr@vn6.1 "jules-vn6.1"

cd jules-vn6.1
export JULES_ROOT=$PWD
echo $JULES_ROOT

# note: JULES already comes with make.cfg. Why does this need to be overwritten?

cp /home/clara/make.cfg ~/MODELS/jules-vn6.1/etc/fcm_make
cd $JULES_ROOT
fcm make -j 2 -f etc/fcm-make/make.cfg --new

rosie go

# For rosie go, you will need to download XQuartz if using Mac or Xming for Windows which you need to run before running Putty

## Add the following to your .bashrc file
echo '#zlib' >> ~/.bashrc
echo 'export PATH=/opt/zlib/:$PATH' >> ~/.bashrc
echo '#openmpi' >> ~/.bashrc
echo 'export PATH=/opt/openmpi/bin:$PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=/opt/openmpi/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
echo '#netcdf' >> ~/.bashrc
echo 'export PATH=/opt/netcdf_par/bin:$PATH' >> ~/.bashrc
echo 'export LD_LIBRARU_PATH=/opt/netcdf_par/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc

### Check installations worked properly
cd ~
echo $SHELL
# This should return "/bin/bash"

xmessage -center hello!
# X11 forwarding should be activated. A little window should appear and say Hello
# If this is not happening, check i) you did use "-X" in the options to log in on your ssh command
# ii) your terminal progra, has been set up to accept X11 Forwarding

cat ~/.bashrc
# Check the .bashrc file and make sure it has been set up correctly eg. a PATH= command should be in there

cylc check-software
# If Cylc is installed, a message saying "Core requirements:OK" should display

rose --version

fcm --version

python --version
# This should be vn>=2.6
