
# CHANGEME: This should be the name of your executable
EXE = slip_problem@EXEEXT@

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS = slip_problem.@OBJEXT@

# CHANGEME: Additional libraries
ADDLIBS =

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS =

##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile Ipopt.   #
##########################################################################

# C Compiler command
CC = @CC@

# C Compiler options
CFLAGS = @CFLAGS@

# additional C Compiler options for linking
CLINKFLAGS = @RPATH_FLAGS@

prefix=@prefix@
exec_prefix=@exec_prefix@

# Include directories
@COIN_HAS_PKGCONFIG_TRUE@INCL = `PKG_CONFIG_PATH=@COIN_PKG_CONFIG_PATH@ @PKG_CONFIG@ --cflags ipopt` $(ADDINCFLAGS)
@COIN_HAS_PKGCONFIG_FALSE@INCL = -I@includedir@/coin-or @IPOPTLIB_CFLAGS@ $(ADDINCFLAGS)

# Linker flags
@COIN_HAS_PKGCONFIG_TRUE@LIBS = `PKG_CONFIG_PATH=@COIN_PKG_CONFIG_PATH@ @PKG_CONFIG@ --libs ipopt` @CXXLIBS@
@COIN_HAS_PKGCONFIG_FALSE@LIBS = -L@libdir@ -lipopt @IPOPTLIB_LFLAGS@ @CXXLIBS@

all: $(EXE)

.SUFFIXES: .c .@OBJEXT@

$(EXE): $(OBJS)
	$(CC) $(CLINKFLAGS) $(CFLAGS) -o $@ $(OBJS) $(ADDLIBS) $(LIBS)

clean:
	rm -rf $(EXE) $(OBJS) ipopt.out

.c.@OBJEXT@:
	$(CC) $(CFLAGS) $(INCL) -c -o $@ $<
