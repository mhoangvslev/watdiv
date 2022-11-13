# WatDiv
Reupload for https://dsg.uwaterloo.ca/watdiv/#download. Their website is incredibly slow.

# Overview
Citation should be formatted as follows:

> [5] G. Aluç, O. Hartig, M. T. Özsu and K. Daudjee. Diversified Stress Testing of RDF Data Management Systems. In Proc. The Semantic Web - ISWC 2014 - 13th International Semantic Web Conference, 2014, pages 197-212.

# Data
Provided that you include a citation to [5], you are free to download and use the WatDiv Data and Query Generator (v0.6) from source code (md5sum=9eac247dfdec044d7fa0141ea3ad361f). The software is supplied "as is" and all use is at your own risk.

Executable files are also provided. Source code and executable files of all versions and changelog can be found [here](https://dsg.uwaterloo.ca/watdiv/changelog.shtml).

The datasets used in the experiments in [5], as well as a billion triples dataset are also available for download:

- [10M Triples](https://dsg.uwaterloo.ca/watdiv/watdiv.10M.tar.bz2) 58,558,746 bytes
- [100M Triples](https://dsg.uwaterloo.ca/watdiv/watdiv.100M.tar.bz2) 629,608,436 bytes
- [1B Triples](https://dsg.uwaterloo.ca/watdiv/watdiv.1000M.tar.bz2) 6,502,656,740 bytes

You may also download the [test workloads](https://dsg.uwaterloo.ca/watdiv/stress-workloads.tar.gz) used in [5].

# Instruction
1. Compiling WatDiv (in C++) is straightforward---the only dependencies are the Boost libraries and the Unix words file (i.e., make sure you have a wordlist package installed under /usr/share/dict/). Once you have installed Boost, simply execute the following commands on UNIX:

```bash
tar xvf watdiv_v05.tar
cd watdiv
setenv BOOST_HOME <BOOST-INSTALLATION-DIRECTORY>
export BOOST_HOME=<BOOST-INSTALLATION-DIRECTORY> (in bash)
make
cd bin/Release
```

2. The last step above is important. To run the data generator, issue the following command:

```bash
./watdiv -d <model-file> <scale-factor>
```
