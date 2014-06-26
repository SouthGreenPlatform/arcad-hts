VERSION="0.1"
INSTALL=${HOME}

all:
	@echo "this package requires no compilation"

check:
	@echo "this package doesn't have any test (yet)"

install:
	mkdir -p ${INSTALL}/bin
	cp sp5_gbs/* ${INSTALL}/bin/

dist:
	mkdir -p arcad-hts-${VERSION}
	cp AUTHORS LICENSE Makefile README TODO arcad-hts-${VERSION}/
	cp -r sp5_gbs/ arcad-hts-${VERSION}/
	tar -czvf arcad-hts-${VERSION}.tar.gz arcad-hts-${VERSION}/
	rm -rf arcad-hts-${VERSION}
