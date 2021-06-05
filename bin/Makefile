CC=g++
INCDIR=../include/
SRC=../src/
MAIN_TEST=../main_test/
CPPFLAG=-I${INCDIR} -std=c++11 -llapacke -lblas

default:run

run:main1 main2 Test1_NoHom Test1_Hom Test1_Hom_Final_RE Test1_NoHom_Final_RE

Test1_NoHom_Final_RE:${MAIN_TEST}Test1_NoHom_Final_RE.cpp ${INCDIR}*.h ${SRC}*.cpp
	g++ ${^}  -g -o ${@} ${CPPFLAG}

Test1_Hom_Final_RE:${MAIN_TEST}Test1_Hom_Final_RE.cpp ${INCDIR}*.h ${SRC}*.cpp
	g++ ${^}  -g -o ${@} ${CPPFLAG}

Test1_Hom:${MAIN_TEST}Test1_Hom.cpp ${INCDIR}*.h ${SRC}*.cpp
	g++ ${^}  -g -o ${@} ${CPPFLAG}

Test1_NoHom:${MAIN_TEST}Test1_NoHom.cpp ${INCDIR}*.h ${SRC}*.cpp
	g++ ${^}  -g -o ${@} ${CPPFLAG}

main1:${MAIN_TEST}main1.cpp ${INCDIR}*.h ${SRC}*.cpp
	g++ ${^}  -g -o ${@} ${CPPFLAG}

main2:${MAIN_TEST}main2.cpp ${INCDIR}*.h ${SRC}*.cpp
	g++ ${^}  -g -o ${@} ${CPPFLAG}

.PHONY:clean

clean:
	-rm main1 main2 -f *.out *.m
