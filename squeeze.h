////////////////////////////////////////////////////////////////
// Compress/decompress functions 
////////////////////////////////////////////////////////////////

// inbuf -- pointer to an array with uncompressed data
// intotal -- length of inbuf in bytes
// outbuf -- destination array (should also be intotal bytes long)
//
// Return values: 0 if memory could not be allocated or inbuf could
// not be compressed. Number of bytes written to outbuf otherwise.
// The chances are you will never see a zero when compressing large
// images.  Just to make sure though, you might want to store an
// extra byte with every image. This byte could be 0 if the image is
// uncompressed, and 1 if it is compressed. If you try to decompress
// an image that is not compress, decompress will not work properly.              
// 
unsigned int compress1(unsigned char *inbuf, unsigned char *outbuf, 
		       long intotal);

// inbuf -- pointer to an array with COMPRESSED data.  If it points
// to anything else, the behavior of expand is unpredictable, segment
// fault being the most likely outcome. 
// outbuf -- destination array. 

void expand(unsigned char* inbuf, unsigned char* outbuf);

// end of squeeze.h
