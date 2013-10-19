CPPFLAGS=-I./ 
SRC=*.c
OUT=plumed.so

default: $(OUT)

links:
	ln -s -f ../common_files/* .

$(OUT): links
	gcc  $(CPPFLAGS) -O3 -DACEMD -fPIC --shared -o $(OUT) *.c

clean:
	rm -f $(OUT) 
	find -type l -delete
	find -name \*~ -delete

veryclean: clean
	rm *~ TAGS
