CC = gcc

CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors 
LDFLAGS = -lm 


SRCS = symnmf.c
OBJS = $(SRCS:.c=.o)

EXEC = symnmf

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJS) $(EXEC)

.PHONY: all clean