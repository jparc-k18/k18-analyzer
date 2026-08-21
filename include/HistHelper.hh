/*
  HistHelper.hh

  2012/1/24
*/

#ifndef HIST_HELPER_H
#define HIST_HELPER_H

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TTree.h>
#include <TString.h>

inline void HB1( int id, const char *title, 
                   int nbinx, double xlow, double xhigh )
{
  new TH1F( Form("h%d",id), title, nbinx, xlow, xhigh );
}

inline void HB2( int id, const char *title,
                   int nbinx, double xlow, double xhigh,
                   int nbiny, double ylow, double yhigh )
{
  new TH2F( Form("h%d",id), title, nbinx, xlow, xhigh,
                   nbiny, ylow, yhigh );
}

inline void HBProf( int id, const char *title, 
                   int nbinx, double xlow, double xhigh , 
                    double ylow, double yhigh)
{
  new TProfile( Form("h%d",id), title, nbinx, xlow, xhigh, ylow, yhigh);
}

inline void HBTree(const char *name, const char *title)
{
  new TTree( name, title);
}

//_________________________________________________________________________
inline void HF1( const char *name, double x, double w=1.0 )
{
  if( TH1 *h = dynamic_cast<TH1 *>(gFile->Get(name)) )
      h->Fill( x, w );
}

inline void HF1( int id, double x, double w=1.0 )
{
  HF1( Form("h%d",id), x, w );
}

//_________________________________________________________________________
inline void HF2( const char *name, double x, double y, double w=1.0 )
{
  if( TH2 *h = dynamic_cast<TH2 *>(gFile->Get(name)) )
    h->Fill( x, y, w );
}

inline void HF2( int id, double x, double y, double w=1.0 )
{
  HF2( Form("h%d",id), x, y, w );
}

//_________________________________________________________________________
inline void HFProf( const char *name, double x, double y, double w=1.0 )
{
  if( TProfile *h = dynamic_cast<TProfile *>(gFile->Get(name)) )
      h->Fill( x, y, w );
}

inline void HFProf( int id, double x, double y, double w=1.0 )
{
  HFProf( Form("h%d",id), x, y, w );
}

#endif
