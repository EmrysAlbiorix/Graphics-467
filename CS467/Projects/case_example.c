

#include <stdio.h>


int testA(int x)
{

  if (x == 1) {
    printf("one\n") ;
  } else if (x == 2) {
    printf("two\n") ;
  } else if (x == 3) {
    printf("three\n") ;
  } else {
    printf("I don't understand\n") ;
  }

}



int testB(int x)
{

  switch (x) {

    case 1 :
      printf("one\n") ;
      break ;

    case 2 :
      printf("two\n") ;
      break ;

    case 3 :
      printf("three\n") ;
      break ;
	
    default :
      printf("I don't understand\n") ;
      break ;
    
  } // end switch

}


int main()
{
  int i ;
  printf("enter an integer : ") ;
  scanf("%d",&i) ;
  testA(i) ;
  testB(i) ;
}
