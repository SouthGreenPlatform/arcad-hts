VERSION=1.0.1# http://semver.org/
INSTALL=${HOME}
PERL_LIB=${INSTALL}/lib
PYTHON_LIB=${INSTALL}/lib
LIBPE=$(shell mkdir -p ${PERL_LIB}; cd ${PERL_LIB}; pwd)
LIBPY=$(shell mkdir -p ${PYTHON_LIB}; cd ${PYTHON_LIB}; pwd)

all:
	@echo "this package requires no compilation"

check:
	@echo "this package doesn't have any test (yet)"

install: 
	mkdir -p ${INSTALL}/bin
	cp scripts/* ${INSTALL}/bin/
	cp sp5_gbs/* ${INSTALL}/bin/
	mkdir -p ${PERL_LIB}/Modules
	cp -r lib/Config ${PERL_LIB}/Modules
	cp -r lib/Files  ${PERL_LIB}/Modules
	cp -r lib/Levenshtein ${PYTHON_LIB}
	cp -r lib/davem_fastq* ${PYTHON_LIB}
	sed -i 's:/NAS/arcad_data/Softs/tags:'$(LIBPE)':' ${INSTALL}/bin/arcad_hts_*.pl
	sed -i 's:version="";:version="'$(VERSION)'";:' ${INSTALL}/bin/arcad_hts_*.pl
	sed -i 's:sys.path".append(""):sys.path.append("'$(LIBPY)'"):' ${INSTALL}/bin/*.py

dist:
	mkdir -p arcad-hts-${VERSION}
	cp AUTHORS LICENSE Makefile README TODO arcad-hts-${VERSION}/
	cp -r scripts arcad-hts-${VERSION}/
	cp -r sp5_gbs/ arcad-hts-${VERSION}/
	sed -i 's:version="";:version="'$(VERSION)'";:' ${INSTALL}/bin/arcad_hts_*.pl
	tar -czvf arcad-hts-${VERSION}.tar.gz arcad-hts-${VERSION}/
	rm -rf arcad-hts-${VERSION}
