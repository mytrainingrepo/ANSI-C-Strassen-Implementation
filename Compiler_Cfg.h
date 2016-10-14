

#ifndef COMPILER_CFG_H_
#define COMPILER_CFG_H_

/* INLINE  define for abstraction of the keyword inline*/
#define INLINE     inline

/* LOCAL_INLINE define for abstraction of the keyword inline in functions with "static" scope.
   Different compilers may require a different sequence of the keywords "static" and "inline"
   if this is supported at all. */
#define LOCAL_INLINE    //static inline

#endif /* COMPILER_CFG_H_ */
