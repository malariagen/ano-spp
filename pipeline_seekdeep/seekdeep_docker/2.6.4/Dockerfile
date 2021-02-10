FROM ubuntu

#Here we follow instructions for de-novo installation from the dev 
# http://seekdeep.brown.edu/installingSeekDeep_Ubuntu.html#brand-new-ubuntu-machine
# Exceptions:
#- do not set timezone
#- do not use sudo as executed as root
#- skip gcc-7 installation from ppa, as apt repos for Ubuntu 18.04 contains these
#- do not add SeekDeep to path, this will be done in Dockerfile

# install compilers
RUN	apt-get update && \
	apt-get dist-upgrade -y && \
	apt-get -y autoremove && \
	apt-get install -y build-essential \
		git \
		make \
		software-properties-common

# clone and install seekdeep
RUN cd /opt/; \
	head -n -1 ~/.profile > tmp.txt; \ 
	mv tmp.txt ~/.profile; \
	echo "test -t 0 && mesg n" >> ~/.profile; \
	git clone https://github.com/bailey-lab/SeekDeep; \
	cd SeekDeep; \
	git checkout tags/2.6.4; \
	./setup.py --libs cmake:3.7.2 --symlinkBin; \
	echo "" >> ~/.profile; \
	echo "#Add SeekDeep bin to your path" >> ~/.profile; \
	echo "export PATH=\"$(pwd)/bin:\$PATH\"" >> ~/.profile; \
	. ~/.profile; \
	./setup.py --addBashCompletion; \
	./install.sh 2; \
	./setup.py --libs muscle:3.8.31 --symlinkBin --overWrite

ENV PATH "/opt/SeekDeep/bin:$PATH"

WORKDIR /data
