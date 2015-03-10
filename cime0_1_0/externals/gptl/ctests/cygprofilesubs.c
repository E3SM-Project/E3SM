#include <stdio.h>

static int junk;
extern void callutil10times ();
extern void callutil100times ();
extern void util ();
extern void A ();
extern void B ();
extern void C ();

void callsubs (int niter)
{
  callutil10times ();
  callutil100times ();
  A();
  util ();
}

void callutil10times ()
{
  int n;
  for (n = 0; n < 10; ++n) 
    util ();
}

void callutil100times ()
{
  int n;
  for (n = 0; n < 100; ++n)
    util ();
}

void util () 
{
  junk = 11;
}

void A ()
{
  B ();
}

void B ()
{
  C ();
}

void C ()
{
  util ();
  callutil10times ();
}
