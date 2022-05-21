CXX=g++

LIBRARIES=-lcmatrix -lrt -lopenblas 
# -lMinuit  -lg2c

# Not generally helpful to include optimisation flags below - optimisation is very time consuming and gains little
# -O3
CXXFLAGS= -I../spirit-1.8.5-miniboost -I../libcmatrix/include -msse3 -mfpmath=sse -DNEED_SEPARATE_COUNT_T -std=c++11 -DNMRSIM_USE_CXX11  -Wall -Wno-sign-compare -DLCM_USE_EXTERNAL -DHAVE_SYS_RESOURCE -DHAVE_STRTOK_R -DUSE_LIBCMATRIX -g  -Wall -Wno-sign-compare
# -DENABLE_DL -DUSE_MINUIT
PNMRSIMHOME= ..

#Used automatically by some makes, otherwise need explicit rule below
COMPILE=$(CXX) $(CXXFLAGS) -c
LINK= $(CXX) $(CXXFLAGS) -L../libcmatrix/lib

%.o: %.cc
	$(COMPILE) -o $@ $<

createdip: createdip.o geometry_utils.o
	$(LINK) -o $@ createdip.o geometry_utils.o $(LIBRARIES)

crys_dip: crys_dip.o geometry_utils.o
	$(LINK) -o $@ crys_dip.o geometry_utils.o $(LIBRARIES)

testtensorrot: testtensorrot.o geometry_utils.o
	$(LINK) -o $@ testtensorrot.o geometry_utils.o $(LIBRARIES)

analyseMQ: analyseMQ.o
	-cd $(PNMRSIMHOME); make stubs
	$(LINK) -o $@ analyseMQ.o $(PNMRSIMHOME)/NMR.o $(PNMRSIMHOME)/Parser.o $(PNMRSIMHOME)/NMR_RFstub.o $(PNMRSIMHOME)/parser_RF.o $(PNMRSIMHOME)/expression_definition.o $(PNMRSIMHOME)/parser_definition.o $(PNMRSIMHOME)/standard_functions.o $(PNMRSIMHOME)/Interaction.o $(PNMRSIMHOME)/fitstub.o $(PNMRSIMHOME)/parser_additional.o $(PNMRSIMHOME)/Loopingstub.o $(PNMRSIMHOME)/NMR_acqstub.o $(PNMRSIMHOME)/Actionstub.o $(PNMRSIMHOME)/common.o $(LIBRARIES)

compspec: compspec.o
	$(LINK) -o $@ compspec.o $(LIBRARIES)

calcaverage: calcaverage.o
	$(LINK) -o $@ calcaverage.o $(LIBRARIES)

quantify: quantify.o
	$(LINK) -o $@ quantify.o $(LIBRARIES)

compfid: compfid.o
	$(LINK) -o $@ compfid.o $(LIBRARIES)

dumpmatlab: dumpmatlab.o
	$(LINK) -o $@ dumpmatlab.o $(LIBRARIES)

Hdensity: Hdensity.o
	$(LINK) -o $@ Hdensity.o $(LIBRARIES)

magres2pNMRsim: magres2pNMRsim.o geometry_utils.o newmagres_parser.o
	$(LINK) -o $@ magres2pNMRsim.o geometry_utils.o newmagres_parser.o $(LIBRARIES)

clean:
	-rm *.o *.exe createdip crys_dip calcaverage compspec quantify compfid dumpmatlab Hdensity magres2pNMRsim testtensorrot testdipolarrot readnewmagres cryswibble cryswibbleold CH2wibble calcav2


