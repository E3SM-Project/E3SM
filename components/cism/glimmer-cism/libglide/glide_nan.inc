#ifndef GLIDE_NAN_INC
#define GLIDE_NAN_INC

#define IS_NAN(x) ((x) /= (x))
#define IS_INF(x) (ISNAN((x)*0))
#define IS_POS_INF(x) (IS_INF(x) .and. ((x) > 0)
#define IS_NEG_INF(x) (IS_INF(x) .and. ((x) < 0)
#endif
