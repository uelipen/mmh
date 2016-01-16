/*
Write gif file.
*/

#ifdef ARCH_IS_CRAY
void WGIF(image,width,height,file)
#else
void wgif_(image,width,height,file)
#endif
int	*width,*height;
char	*image,*file;
{
	wgif(image,*width,*height,file);
}

#include <stdio.h>
#define BITS_PER_PIXEL	8

int wgif(image,width,height,file)

int	width,height;
char	*image,*file;

{
	int	i,codesize,lenf=80;
	char	fbuf[255];
	FILE	*fp;
	extern void fputshort(),compress();
#include "cmap.h"

	while (--lenf > 0 && (file[lenf]==' ' || file[lenf]==0));
	for (i=0; i<=lenf; i++) fbuf[i] = file[i];
	fbuf[lenf+1] = '\0';
	fp = fopen(fbuf,"w");
	if (!fp) {
		fprintf(stderr,"Cannot open file %s\n",file);
		return(-1);
	}
	fwrite("GIF87a",1,6,fp);/* GIF signature. */

/* Screen Descriptor section. */

	fputshort(width,fp);	/* Screen width. */
	fputshort(height,fp);	/* Screen height. */
	i = ( 0x80 |  (BITS_PER_PIXEL-1) | ((BITS_PER_PIXEL-1)*16));
	fputc(i,fp);		/* Color descriptor byte. */
	fputc(0,fp);		/* Background color byte. */
	fputc(0,fp);		/* Zero byte to end screen descriptor. */

	for (i=0;i<768;i++)
		fputc(cmap[i],fp);	/* Color map. */

/* Image Descriptor section. */
	fputc(0x2c,fp);		/* Image separator character. */
	i = 0;
	fputshort(i,fp);	/* Image left. */
	fputshort(i,fp);	/* Image top. */
	fputshort(width,fp);	/* Image width. */
	fputshort(height,fp);	/* Image height. */
	fputc(0,fp);		/* Color descriptor byte: 0x00 for global color map, sequential ordering. */

/* Raster data section. */
	codesize = BITS_PER_PIXEL;
	fputc(codesize,fp);	/* Initial code size */
	compress(codesize+1,fp,image,width*height);
	fputc(0,fp);		/* End of data stream */

	fputc(0x3b,fp);		/* GIF terminator byte. */

	if (ferror(fp)) {
		fprintf(stderr,"There was an error writing the GIF file\n");
		return(-1);
	}
	fclose(fp);
	return(0);

}

void fputshort(word, fp)

int	word;
FILE	*fp;

{
  /* writes a 16-bit integer in GIF order (LSB first) */
  fputc((word&0xff),fp);
  fputc(((word>>8)&0xff),fp);
}
