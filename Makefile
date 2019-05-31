#CDEBUGFLAGS = -g
CDEBUGFLAGS = -O2
CFLAGS = ${CDEBUGFLAGS} -I ${HOME}/work/obj
CC = gcc

#----------------------------------------------------------*

NM_HDRS_ALL = nr.h
NM_CSRS_ALL = nr.c
NM_OBJS_ALL = nr.o

rlf: rlf.o rlfn.o ${NM_OBJS_ALL}
	${CC} ${CFLAGS} rlf.o rlfn.o ${NM_OBJS_ALL} -lm -o rlf.ex

rlf.o:  rlf.c rlf.h rlfa.h
	${CC} ${CFLAGS} -c rlf.c
rlfn.o: rlfn.c rlfn.h rlfa.h
	${CC} ${CFLAGS} -c rlfn.c
nr.o:  nr.c nr.h
	${CC} ${CFLAGS} -c nr.c

clean:
	/bin/rm -f *.o

.KEEP_STATE:
 
