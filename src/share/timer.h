#ifdef TIMING
#define TIMER_START(c)  c=CLOCK()
#define TIMER_DETAIL_START(t,d,c)  if(d<t%detail) c=CLOCK()
#define TIMER_UPDATE(t,i,e,s) t%time_buf(i)=t%time_buf(i)+e-s
#define TIMER_SYNC(a,b)  a=b	
#define TIMER_DETAIL(i,t) (i<t%detail)
#else
#define TIMER_START(c)  c=0
#define TIMER_DETAIL_START(t,d,c) !   
#define TIMER_UPDATE(t,i,e,s) s=e
#define TIMER_SYNC(a,b) a=b
#define TIMER_DETAIL(i,t) (.false.)
#endif
