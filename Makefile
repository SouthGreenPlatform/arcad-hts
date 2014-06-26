VERSION="0.1"

all:
	@echo "this package requires no compilation"

check:
	@echo "this package doesn't have any test (yet)"

dist:
	mkdir -p arcad-hts-${VERSION}
	cp sp5_gbs/* arcad-hts-${VERSION}
	tar -czvf arcad-hts-${VERSION}.tar.gz arcad-hts-${VERSION}/
	rm -rf arcad-hts-${VERSION}
