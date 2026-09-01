# plot() draws. Without a null device every drawing test would open the
# default one and leave an Rplots.pdf behind in the package sources, which
# then gets swept into the source tarball. pdf(NULL) writes nothing and is
# closed when the session ends, so no teardown is needed.
grDevices::pdf(NULL)
