#include <stdio.h>
#include <stdlib.h>
#include "squeeze.h"

#define BITS 14
                   
#define HASHING_SHIFT BITS-8      
#define MAX_VALUE (1 << BITS) - 1 
#define MAX_CODE MAX_VALUE - 1    
                                  
#if BITS == 14
#   define TABLE_SIZE 18041        /* The string table size needs to be a */
#endif                            /* prime number that is somewhat larger*/
#if BITS == 13                    /* than 2**BITS.                       */
#   define TABLE_SIZE 9029
#endif
#if BITS <= 12
#   define TABLE_SIZE 5021
#endif

// These are global pointers to input and output arrays

unsigned int inbufp;
unsigned int outbufp;

int *code_value;                  /* This is the code value array        */
unsigned int *prefix_code;        /* This array holds the prefix codes   */

unsigned char *append_character;  /* This array holds the appended chars */
unsigned char decode_stack[4000]; /* This array holds the decoded string */


unsigned int  find_match(int hash_prefix, unsigned int hash_character);
int     output_code(unsigned char* outbuf, long intotal, unsigned int code);
unsigned int  input_code(unsigned char* inbuf);
unsigned char *decode_string(unsigned char *buffer, unsigned int code);

void fcompress_(unsigned char *inbuf, unsigned char *outbuf, 
		int *outtotal, int *intotal)
{
	long lintotal;
	lintotal=*intotal;
	*outtotal=compress1(inbuf,outbuf,lintotal);
}

unsigned int compress1(unsigned char *inbuf, unsigned char *outbuf, 
		       long intotal)
{
unsigned int next_code;
unsigned int character;
unsigned int string_code;
unsigned int index;
int i;

  inbufp = 0;
  outbufp = 0;
  intotal--;

  code_value= (int *) malloc(TABLE_SIZE*sizeof(int));
  prefix_code=  (unsigned int *) malloc(TABLE_SIZE*sizeof(int));
  append_character= (unsigned char *)  malloc(TABLE_SIZE*sizeof(char));

  if (code_value==NULL || prefix_code==NULL || append_character==NULL)
  {
     printf("Compress: not enough memory\n");
     return 0;  
  }

  next_code=256;              /* Next code is the next available string code*/
  for (i=0;i<TABLE_SIZE;i++)  /* Clear out the string table before starting */
    code_value[i]=-1;

  i=0;
  string_code=inbuf[inbufp++];    /* Get the first code      */

  /*
  ** This is the main loop where it all happens.  This loop runs util all of
  ** the input has been exhausted.  Note that it stops adding codes to the
  ** table after all of the possible codes have been defined.
  */

  while (inbufp <= intotal)
    {
    character=inbuf[inbufp++];

    index=find_match(string_code,character);  /* See if the string is in */
    if (code_value[index] != -1)              /* the table.  If it is,   */
      string_code=code_value[index];          /* get the code value.  If */
    else                                      /* the string is not in the*/
      {                                       /* table, try to add it.   */
      if (next_code <= MAX_CODE)
	{
	  code_value[index]=next_code++;
	  prefix_code[index]=string_code;
	  append_character[index]=character;
	}
      output_code(outbuf, intotal, string_code);
      string_code=character;                                
      }                                                     
    }                                                       
 
  /*
  ** End of the main loop.
  */

  output_code(outbuf, intotal, string_code); /* Output the last code                */
  output_code(outbuf, intotal, MAX_VALUE);   /* Output the end of buffer code       */
  i = output_code(outbuf, intotal, 0);       /* This code flushes the output buffer */

  free(code_value);
  free(prefix_code);
  free(append_character);

  if (i == 0)
     return 0;
  else 
     return outbufp;
}

/*
** This is the hashing routine.  It tries to find a match for the prefix+char
** string in the string table.  If it finds it, the index is returned.  If
** the string is not found, the first available index in the string table is
** returned instead.
*/

unsigned int find_match(int hash_prefix, unsigned int hash_character)
{
int index;
int offset;

  index = (hash_character << HASHING_SHIFT) ^ hash_prefix;
  if (index == 0)
    offset = 1;
  else
    offset = TABLE_SIZE - index;
  while (1)
    {
    if (code_value[index] == -1)
      return(index);
    if (prefix_code[index] == hash_prefix &&
	append_character[index] == hash_character)
      return(index);
    index -= offset;
    if (index < 0)
      index += TABLE_SIZE;
    }
}

/*
**  This is the expansion routine.  It takes an LZW format file, and expands
**  it to an output file.  The code here should be a fairly close match to
**  the algorithm in the accompanying article.
*/

void expand(unsigned char* inbuf, unsigned char* outbuf)
{
unsigned int next_code;
unsigned int new_code;
unsigned int old_code;
int character;
int counter;
unsigned char *string;

  prefix_code= (void *) malloc(TABLE_SIZE*sizeof(int));
  append_character= (void *) malloc(TABLE_SIZE*sizeof(char));
  
  if (prefix_code==NULL || append_character==NULL)
    {
     printf("Decompress: not enough memory\n");
     return;  
    }

  next_code=256;           /* This is the next available code to define */
  counter=0;               /* Counter is used as a pacifier.            */

  inbufp = 0;
  outbufp = 0;

  old_code=input_code(inbuf);         /* Read in the first code, initialize the */
  character=old_code;                 /* character variable, and send the first */
  outbuf[outbufp++] = old_code;       /* code to the output file                */
  
  /*
  **  This is the main expansion loop.  It reads in characters from the LZW file
  **  until it sees the special code used to inidicate the end of the data.
  */

  while ((new_code=input_code(inbuf)) != (MAX_VALUE))
    {
      /*
      ** This code checks for the special STRING+CHARACTER+STRING+CHARACTER+STRING
      ** case which generates an undefined code.  It handles it by decoding
      ** the last code, and adding a single character to the end of the decode string.
      */
    if (new_code>=next_code)
      {
      *decode_stack = character;
      string=decode_string(decode_stack+1,old_code);
      }
    /*
    ** Otherwise we do a straight decode of the new code.
    */
    else
      string=decode_string(decode_stack,new_code);
    /*
    ** Now we output the decoded string in reverse order.
    */
    character=*string;
    while (string >= decode_stack)
      outbuf[outbufp++] = *(string--);
    /*
    ** Finally, if possible, add a new code to the string table.
    */
    if (next_code <= MAX_CODE)
      {
      prefix_code[next_code]=old_code;
      append_character[next_code]=character;
      next_code++;
      }
    old_code=new_code;
    }

  free(prefix_code);
  free(append_character);
}

/*
** This routine simply decodes a string from the string table, storing
** it in a buffer.  The buffer can then be output in reverse order by
** the expansion program.
*/

unsigned char *decode_string(unsigned char *buffer,unsigned int code)
{
int i;

  i=0;
  while (code > 255)
    {
    *(buffer++) = append_character[code];
    code=prefix_code[code];
    if (i++>=4094)
      {
      printf("Fatal error during code expansion.\n");
      exit(1);
      }
    }
  *buffer=code;
  return(buffer);
}

/*
** The following two routines are used to output variable length
** codes.  They are written strictly for clarity, and are not
** particularyl efficient.
*/

unsigned int input_code(unsigned char* inbuf)
{
unsigned int return_value;
static int input_bit_count=0;
static unsigned int input_bit_buffer=0;

  while (input_bit_count <= 24)
    {
    input_bit_buffer |=
      (unsigned int) inbuf[inbufp++] << (24-input_bit_count);
    input_bit_count += 8;
    }
  return_value=input_bit_buffer >> (32-BITS);
  input_bit_buffer <<= BITS;
  input_bit_count -= BITS;
  return(return_value);
}

int output_code(unsigned char* outbuf, long intotal, unsigned int code)
{
static int output_bit_count=0;
static unsigned int output_bit_buffer=0;

  if (outbufp >= intotal)
    return 0;  

  output_bit_buffer |= (unsigned int) code << (32-BITS-output_bit_count);
  output_bit_count += BITS;
  while (output_bit_count >= 8) 
    {
    outbuf[outbufp++] = output_bit_buffer >> 24;
    output_bit_buffer <<= 8;
    output_bit_count -= 8;
    }
  return 1;
}






