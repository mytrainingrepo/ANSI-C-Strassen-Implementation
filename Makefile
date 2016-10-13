
SRCS1= start.c
SRCS2= memleak.c
CFLAGS= -g -O0 
TARGET1= start
TARGET2= memleak

$(TARGET1): $(SRCS1)
	$(CC) $(CFLAGS) -o $(TARGET1).exe $(SRCS1)

$(TARGET2): $(SRCS2)
	$(CC) $(CFLAGS) -o $(TARGET2).exe $(SRCS2)

clean:
	rm -f *.o
	rm -f *.exe
	rm -f *.a
	rm -f *.so.*
	rm -f *.out
