VERSION	       = v0.3.2
override LIBS += -lz

ifneq ($(wildcard $(ROOTSYS)/lib/root),)
        ROOTLIBS = -L$(ROOTSYS)/lib/root -lCore -lRIO -lHist -lGraf -lGpad -lTree -lMathCore
else
        ROOTLIBS = -L$(ROOTSYS)/lib      -lCore -lRIO -lHist -lGraf -lGpad -lTree -lMathCore
endif

ifneq ($(wildcard $(ROOTSYS)/include/root),)
        INC = -I$(ROOTSYS)/include/root -I$(SAMDIR)
else
        INC = -I$(ROOTSYS)/include      -I$(SAMDIR)
endif

SAMDIR = samtools
SAMLIB = $(SAMDIR)/libbam.a
HTSDIR = $(wildcard $(SAMDIR)/htslib-*)
ifneq ($(HTSDIR),)
        SAMLIB += $(HTSDIR)/libhts.a
        INC    += -I$(HTSDIR)
endif

ifeq ($(OMP),no)
        $(info Compiling with NO parallel support)
else
        OMPFLAGS = -fopenmp
        $(info Compiling with parallel (OpenMP) support)
endif

ifneq ($(YEPPPLIBDIR),)
        override LIBS += -L$(YEPPPLIBDIR) -lyeppp
endif

ifneq ($(YEPPPINCLUDEDIR),)
        INC += -I$(YEPPPINCLUDEDIR) -DUSE_YEPPP
endif

CXX    = g++ -O3 -std=c++11 -DCNVNATOR_VERSION=\"$(VERSION)\" $(OMPFLAGS)

OBJDIR = obj
OBJS   = $(OBJDIR)/cnvnator.o  \
	 $(OBJDIR)/EXOnator.o  \
	 $(OBJDIR)/HisMaker.o  \
	 $(OBJDIR)/AliParser.o \
	 $(OBJDIR)/Genotyper.o \
	 $(OBJDIR)/Interval.o  \
	 $(OBJDIR)/Genome.o

DISTRIBUTION = $(PWD)/CNVnator_$(VERSION).zip
TMPDIR	     =  /tmp
CNVDIR	     = CNVnator_$(VERSION)
MAINDIR	     = $(TMPDIR)/$(CNVDIR)
SRCDIR	     = $(MAINDIR)/src

all: cnvnator

cnvnator: $(OBJS)
	$(CXX) -o $@ $(OBJS) $(SAMLIB) $(LIBS) $(ROOTLIBS)

$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(OBJDIR)
	$(CXX) $(INC) -c $< -o $@

clean:
	rm -fr $(OBJDIR) cnvnator

distribution: clean all
	@echo Creating directory ...
	@rm -rf $(MAINDIR)
	@rm -f  $(DISTRIBUTION)
	@mkdir  $(MAINDIR)
	@mkdir  $(SRCDIR)
	@echo Copying files ...
	@cp *.hh *.cpp  $(SRCDIR)
	@cp Makefile    $(SRCDIR)
	@cp -r samtools $(SRCDIR)
	@rm -f $(SRCDIR)/samtools/samtools
	@rm -f $(SRCDIR)/samtools/*.o
	@rm -f $(SRCDIR)/samtools/*/*.o
	@rm -f $(SRCDIR)/samtools/*/*/*.o
	@rm -f $(SRCDIR)/samtools/*.a
	@rm -f $(SRCDIR)/samtools/*/*.a
	@rm -f $(SRCDIR)/samtools/*/*/*.a
	@rm -fr $(SRCDIR)/samtools/.svn
	@rm -fr $(SRCDIR)/samtools/*/.svn
	@rm -fr $(SRCDIR)/samtools/*/*/.svn
	@cp README          $(MAINDIR)
	@cp CITATION        $(MAINDIR)
	@cp license.rtf     $(MAINDIR)
	@cp cnvnator2VCF.pl $(MAINDIR)
	@echo Zipping ...
	@ln -s $(MAINDIR)
	@zip -qr $(DISTRIBUTION) $(CNVDIR)
	@rm $(CNVDIR)
