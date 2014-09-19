/*
 * Pstools.h
 */

#ifndef PSTOOLS_H_
#define PSTOOLS_H_


#include "TH1F.h"
#include "TPad.h"
#include "TString.h"

class PSTools {
 public:

	static void Pcolor(int icol = 1);
	static void white()   {Pcolor(0);}
	static void black()   {Pcolor(1);}
	static void red()     {Pcolor(2);}
	static void green()   {Pcolor(3);}
	static void blue()    {Pcolor(4);}
	static void yellow()  {Pcolor(5);}
	static void magenta() {Pcolor(6);}
	static void turqouis(){Pcolor(7);}
	static void green2()  {Pcolor(8);}
	static void blue2()   {Pcolor(9);}

	static void coutf(Int_t w,TString s);
	static void coutf(Int_t w,Int_t i);
	static void coutf(int w,std::string s);
	static void coutf(int w,double x);
	static void coutf(int w,int p,double x);

};

#endif /* PSTOOLS_H_ */
