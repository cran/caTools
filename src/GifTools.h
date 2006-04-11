#ifndef GIF_TOOLS_H
#define GIF_TOOLS_H

extern "C" {
  #include <R.h>
  #include <Rinternals.h> 
  #define print Rprintf
  #define Error error
  typedef unsigned char uchar;

  int imreadGif(char* filename, int nImage, bool verbose,
              uchar** data, int &nRow, int &nCol, int &nBand,
              int ColorMap[255], int &Transparent, char** Comment);
              
  int imwriteGif(char* filename, const uchar* data, int nRow, int nCol,
                  int nBand, int nColor, const int *ColorMap,  bool interlace, 
                 int transparent, int DalayTime, char* comment);
}
#endif
              
