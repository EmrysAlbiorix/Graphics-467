#include <stdio.h>
#include <stdlib.h>
#include <string.h>



// convert an xwd file that was constructed using the
// linux "convert" utility 
// to an xwd file that is consistent with jeff's usual
// tool set





//static void fputintB (int x, FILE *fp)
static void fputint (int x, FILE *fp)
{
   char *p ;
   int c[4],i ;

   p = (char *) &x ;

   for (i = 0 ; i < 4 ; i++) {
        c[i] = *p ;
        p++ ;
   }

   for (i = 3 ; i >= 0 ; i--) {
        fputc(c[i], fp) ;
   }

}




//static int fgetintB (FILE *fp)
static int fgetint(FILE *fp)
{
   char *p ;
   int c[4],i ;
   int x ;
   char q ;

   for (i = 3 ; i >= 0 ; i--) {
        c[i] = fgetc(fp) ;
   }



   p = (char *) &x ;

   for (i = 0 ; i < 4 ; i++) {
        q = c[i] ;
        *p = q ;;
        p++ ;
   }

   return x ;
}








void convertthem(FILE *inf, FILE *outf)
{
  // assumes full color but uses full 32 bits per pixel

 int header_size, file_version, format, depth ;
 int xoffset, byte_order, bitmap_unit, bitmap_order, bitmap_pad ;
 int bits_per_pixel, bytes_per_line, visual_class ;
 int redmask, greenmask, bluemask ;
 int bitsperrgb, numentriesincolormap, numcolorstructures ;
 int windowwidth, windowheight, windowx, windowy, windowborderwidth ;
 int numbytestowrite ;

 int width,height ;

 int c, numbytestoread ;
 char qq ;
 int ii ;
 char *sourcedata, *resdata ;
 int ss,rr,jj ;
 char aa,bb,cc,dd ;


 header_size = fgetint(inf) ;
 if (header_size != 104) {
      printf("header_size changed from %d to 104\n",header_size) ;
 }
 header_size = 104 ;
 fputint(header_size, outf) ;


 file_version = fgetint(inf) ;
 if (file_version != 7) {
      printf("file_version changed from %d to 7\n",file_version) ;
 }
 file_version = 7 ;
 fputint(file_version, outf) ;


 format = fgetint(inf) ;
 if (format != 2) {
      printf("format changed from %d to 2\n",format) ;
 }
 format = 2 ;
 fputint(format, outf) ;


 depth = fgetint(inf) ;
 if (depth != 24) {
      printf("depth changed from %d to 24\n",depth) ;
 }
 depth = 24 ;
 fputint(depth, outf) ;


 width =  fgetint(inf) ;
printf("width = %d\n",width) ;
 fputint(width, outf) ;


 height =  fgetint(inf) ;
printf("height = %d\n",height) ; 
 fputint(height, outf) ;


 xoffset = fgetint(inf) ;
 if (xoffset != 0) {
      printf("xoffset changed from %d to 0\n",xoffset) ;
 }
 xoffset = 0 ;
 fputint(xoffset, outf) ;


 byte_order = fgetint(inf) ;
 if (byte_order != 0) {
      printf("byte_order changed from %d to 0\n",byte_order) ;
 }
 byte_order = 0 ;
 fputint(byte_order, outf) ;


 bitmap_unit = fgetint(inf) ;
 if (bitmap_unit != 32) {
      printf("bitmap_unit changed from %d to 32\n",bitmap_unit) ;
 }
 bitmap_unit = 32 ;
 fputint(bitmap_unit, outf) ;


 bitmap_order = fgetint(inf) ;
 if (bitmap_order != 0) {
      printf("bitmap_order changed from %d to 0\n",bitmap_order) ;
 }
 bitmap_order = 0 ;
 fputint(bitmap_order, outf) ;


 bitmap_pad = fgetint(inf) ;
 if (bitmap_pad != 32) {
      printf("bitmap_pad changed from %d to 32\n",bitmap_pad) ;
 }
 bitmap_pad = 32 ;
 fputint(bitmap_pad, outf) ;


 bits_per_pixel = fgetint(inf) ;
 if (bits_per_pixel != 32) {
      printf("bits_per_pixel changed from %d to 32\n",bits_per_pixel) ;
 }
 int new_bits_per_pixel = 32 ;
 fputint(new_bits_per_pixel, outf) ;




 
 bytes_per_line = fgetint(inf) ;
 if (bytes_per_line != width*4) {
      printf("bytes_per_line changed from %d to %d\n",bytes_per_line, width*4) ;
 }
 int new_bytes_per_line = width*4 ;
 fputint(new_bytes_per_line, outf) ;

 int pad_bytes_per_input_line = 0 ; 
 if (bits_per_pixel == 24) {
   pad_bytes_per_input_line = bytes_per_line - width*3 ; // for later use
 }


 



 
 

 visual_class = fgetint(inf) ;
 if (visual_class != 5) {
      printf("visual_class changed from %d to 5\n",visual_class) ;
 }
 visual_class = 5 ;
 fputint(visual_class, outf) ;


 redmask = fgetint(inf) ;
 if (redmask != 0x00ff0000) {
      printf("redmask changed from %d to 0x00ff0000\n",redmask) ;
 }
 redmask = 0x00ff0000 ;
 fputint(redmask, outf) ;


 greenmask = fgetint(inf) ;
 if (greenmask != 0x0000ff00) {
      printf("greenmask changed from %d to 0x0000ff00\n",greenmask) ;
 }
 greenmask = 0x0000ff00 ;
 fputint(greenmask, outf) ;


 bluemask = fgetint(inf) ;
 if (bluemask != 0x000000ff) {
      printf("bluemask changed from %d to 0x000000ff\n",bluemask) ;
 }
 bluemask = 0x000000ff ;
 fputint(bluemask, outf) ;


 bitsperrgb = fgetint(inf) ;
 if (bitsperrgb != 24) {
      printf("bitsperrgb changed from %d to 24\n",bitsperrgb) ;
 }
 bitsperrgb = 24 ;
 fputint(bitsperrgb, outf) ;


 numentriesincolormap = fgetint(inf) ;
 if (numentriesincolormap != 0) {
      printf("numentriesincolormap changed from %d to 0\n",numentriesincolormap) ;
 }
 fputint(0, outf) ; // there will be no color map in the new form


 numcolorstructures = fgetint(inf) ;
 if (numcolorstructures != 0) {
      printf("numcolorstructures changed from %d to 0\n",numcolorstructures) ;
 }
 numcolorstructures = 0 ;
 fputint(numcolorstructures, outf) ;


 windowwidth = fgetint(inf) ;
 if (windowwidth != width) {
      printf("windowwidth changed from %d to %d\n",windowwidth,width) ;
 }
 windowwidth = width ;
 fputint(windowwidth, outf) ;


 windowheight = fgetint(inf) ;
 if (windowheight != height) {
      printf("windowheight changed from %d to %d\n",windowheight,height) ;
 }
 windowheight = height ;
 fputint(windowheight, outf) ;


 windowx = fgetint(inf) ;
 if (windowx != 0) {
      printf("windowx changed from %d to 0\n",windowx) ;
 }
 windowx = 0 ;
 fputint(windowx, outf) ;


 windowy = fgetint(inf) ;
 if (windowy != 0) {
      printf("windowy changed from %d to 0\n",windowy) ;
 }
 windowy = 0 ;
 fputint(windowy, outf) ;


 windowborderwidth = fgetint(inf) ;
 if (windowborderwidth != 0) {
      printf("windowborderwidth changed from %d to 0\n",windowborderwidth) ;
 }
 windowborderwidth = 0 ;
 fputint(windowborderwidth, outf) ;





 
 char n[2000] ;
 int k = 0 ;

 // read null terminated window name
 while ((n[k++] = c = fgetc(inf)) != '\0') ;
 n[k] = 0 ;
 printf("window name = .%s.\n",n) ;

 
 // for the jeff version, replace the name with 4 null chars
 c = '\0' ;
 fputc(c,outf) ; fputc(c,outf) ; fputc(c,outf) ; fputc(c,outf) ;



 // now copy over the rest of it


 sourcedata = (char *)malloc(width*3) ;
 if (sourcedata == NULL) {
      printf("ERROR: can't malloc space needed\n") ;
      printf("Program terminating\n\n") ;
      exit(1) ;
 }



 resdata = (char *)malloc(width*4) ;
 if (resdata == NULL) {
      printf("ERROR: can't malloc space needed\n") ;
      printf("Program terminating\n\n") ;
      exit(1) ;
 }








 for (ii = 0 ; ii < height ; ii++) {

    ss = 0 ; rr = 0 ;

    for (jj = 0 ; jj < width ; jj++) {

        aa = fgetc(inf) ;
        bb = fgetc(inf) ;
        cc = fgetc(inf) ;

	fputc(cc,outf) ;
	fputc(bb,outf) ;
	fputc(aa,outf) ;
	fputc( 0,outf) ;	
    }

    // read and discard the padding at the end of the line
    for (k = 0 ; k < pad_bytes_per_input_line  ; k++) {
      fgetc(inf) ;
    }
    
    
    

 }


}





/*

 for (ii = 0 ; ii < height ; ii++) {


    fread (sourcedata, width*3, 1, inf) ;

    ss = 0 ; rr = 0 ;

    for (jj = 0 ; jj < width ; jj++) {

        aa = sourcedata[ss++] ;
        bb = sourcedata[ss++] ;
        cc = sourcedata[ss++] ;

        resdata[rr++] = cc ;
        resdata[rr++] = bb ; 
        resdata[rr++] = aa ; 
        resdata[rr++] = 0 ;
    }
    
    

    fwrite (resdata, width*4, 1, outf) ;

*/






void main(int argc, char **argv)
{
  FILE *inf, *outf ;

  if (argc != 3) {

     printf("usage : pgmname  name_of_XWD_generated_xwd_file  name_of_desired_JEFF_xwd_file\n") ;
     exit(1) ;

  }

 

  inf = fopen(argv[1],"r") ;
  if (inf == NULL) {
      printf("can't open file %s\n",argv[1]) ;
      exit(1) ;
  }  



  outf = fopen(argv[2],"w") ;
  if (outf == NULL) {
      printf("can't open file %s\n",argv[2]) ;
      exit(1) ;
  }  


  convertthem(inf,outf) ;


}



