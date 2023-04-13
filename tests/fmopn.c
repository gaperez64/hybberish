#include <stdio.h>
#include <string.h>


static char buffer[] = "foobar";


int
main (void)
{
    int ch;
    FILE *stream;


    stream = fmemopen(buffer, strlen (buffer), "r");
    if (stream == NULL)
        fprintf(stderr, "got null!\n");


    while ((ch = fgetc(stream)) != EOF)
        fprintf(stderr, "Got %c\n", ch);


    fclose(stream);
    return (0);
}
