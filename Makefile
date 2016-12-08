main:main.o
	gcc -o $@ $^
.c.o:
	gcc -c $<

clean:
	rm main *.o
