class Base
{
 public:
  Base ();
  ~Base ();
};

Base::Base ()
{
}

Base::~Base ()
{
}

class X: Base
{
 public:
  X () 
  {
  }
  ~X() 
  {
  }
  void func (int x)
  {
  }
  void func (double x)
  {
  }
};

class Y: Base
{
 public:
  Y ();
  ~Y();
  void func (int);
  void func (double);
};

Y::Y ()
{
}

Y::~Y()
{
}

void Y::func (int x)
{
}

void Y::func (double x)
{
}
