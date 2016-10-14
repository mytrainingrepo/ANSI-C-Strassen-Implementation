
SRCS= Test.c MatrixOps.c
CFLAGS= -O4 --static
TARGET= Test

$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -o $(TARGET).exe $(SRCS)

clean:
	rm -f $(TARGET).exe
