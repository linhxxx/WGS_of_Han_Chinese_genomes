# customization;  change installDir to suit your needs

installDir = ~/${MACHTYPE}/bin

# default targets

default: install

# installation;  change installDir above to suit your needs

programs = lav_sort gfa_to_r_dotplot

install:
	for p in $(programs); do ../tools/install_py.py ${installDir} $${p} || exit 1; done

