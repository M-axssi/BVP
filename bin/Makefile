CC=g++
INCDIR=../include/
SRC=../src/
MAIN_TEST=../main_test/
CPPFLAG=-I${INCDIR} -std=c++11 -llapacke -lblas

default:run

run:main1 main2

main1:${MAIN_TEST}main1.cpp ${INCDIR}*.h ${SRC}*.cpp
	g++ ${^}  -g -o ${@} ${CPPFLAG}

main2:${MAIN_TEST}main2.cpp ${INCDIR}*.h ${SRC}*.cpp
	g++ ${^}  -g -o ${@} ${CPPFLAG}

.PHONY:clean

clean:
	-rm main1 main2 -f *.out *.m
