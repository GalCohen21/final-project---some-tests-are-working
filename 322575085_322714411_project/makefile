CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
TARGET = symnmf
OBJS = symnmf.o symnmf_main.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) -lm

symnmf.o: symnmf.c symnmf.h
	$(CC) $(CFLAGS) -c symnmf.c

symnmf_main.o: symnmf_main.c symnmf.h
	$(CC) $(CFLAGS) -c symnmf_main.c

clean:
	rm -f $(OBJS) $(TARGET)