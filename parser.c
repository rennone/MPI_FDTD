#include "parser.h"

bool parser_nextLine(FILE *fp, char buf[])
{
  while( fgets(buf, 1024, fp) != NULL)
  {
    //空行はとばす
    if(buf[0] == '#' || buf[0] == '\0' || buf[0] == '\n')
      continue;
    
    //#より後ろはカット
    char* p = strstr(buf, "#"); //#の位置を探す
    if( p != NULL){
      strncpy(buf, buf, p-buf+1); //#の手前までをとってくる
    }
    return true;
  }
  
  return false;
}
