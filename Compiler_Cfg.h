

#ifndef COMPILER_CFG_H_
#define COMPILER_CFG_H_


/* INLINE  define for abstraction of the keyword inline*/
#define INLINE     //inline

/* LOCAL_INLINE define for abstraction of the keyword inline in functions with "static" scope.
   Different compilers may require a different sequence of the keywords "static" and "inline"
   if this is supported at all. */
#define LOCAL_INLINE    //static inline

/* AUTOMATIC used for the declaration of local pointers */
#define AUTOMATIC


/* FUNC_P2CONST macro for declaration and definition of functions returning a pointer to a constant, that ensures
     correct syntax of function declarations.
   rettype     return type of the function
   ptrclass    defines the classification of the pointer�s distance
   memclass    classification of the function itself*/
# define FUNC_P2CONST(rettype, ptrclass, memclass) const rettype ptrclass * memclass

/* FUNC_P2VAR macro for the declaration and definition of functions returning a pointer to a variable, that ensures
     correct syntax of function declarations
   rettype     return type of the function
   ptrclass    defines the classification of the pointer�s distance
   memclass    classification of the function itself*/
# define FUNC_P2VAR(rettype, ptrclass, memclass) rettype ptrclass * memclass

/* P2VAR macro for the declaration and definition of pointers in RAM, pointing to variables
   ptrtype     type of the referenced variable memclass
   memclass    classification of the pointer�s variable itself
   ptrclass    defines the classification of the pointer�s distance
 */

/* FUNC macro for the declaration and definition of functions, that ensures correct syntax of function declarations
   rettype     return type of the function
   memclass    classification of the function itself*/
# define FUNC(rettype, memclass) rettype memclass




/* Type definition of const pointers to functions
   rettype     return type of the function
   ptrclass    defines the classification of the pointer's distance
               (not used on 32bit platforms)
   fctname     function name respectivly name of the defined type
 */
#define CONSTP2FUNC(rettype, ptrclass, fctname)  rettype (*const fctname)


# define P2VAR(ptrtype, memclass, ptrclass) memclass ptrtype ptrclass *

/* P2CONST macro for the declaration and definition of pointers in RAM pointing to constants
   ptrtype     type of the referenced data
   memclass    classification of the pointer's variable itself
   ptrclass    defines the classification of the pointer's distance
 */
# define P2CONST(ptrtype, memclass, ptrclass) const memclass ptrtype ptrclass *

/* CONSTP2VAR macro for the declaration and definition of constant pointers accessing variables
   ptrtype     type of the referenced data
   memclass    classification of the pointer's variable itself
   ptrclass    defines the classification of the pointer's distance
 */
# define CONSTP2VAR(ptrtype, memclass, ptrclass) memclass ptrtype ptrclass * const

/* CONSTP2CONST macro for the declaration and definition of constant pointers accessing constants
   ptrtype     type of the referenced data
   memclass    classification of the pointer's variable itself
   ptrclass    defines the classification of the pointer's distance
 */
# define CONSTP2CONST(ptrtype, memclass, ptrclass) const memclass ptrtype ptrclass * const

/* P2FUNC macro for the type definition of pointers to functions
   rettype     return type of the function
   ptrclass    defines the classification of the pointer's distance
   fctname     function name respectivly name of the defined type
 */
# define P2FUNC(rettype, ptrclass, fctname) rettype (ptrclass * fctname)

/* CONST macro for the declaration and definition of constants
   type        type of the constant
   memclass    classification of the constant itself
 */
//# define CONST(type, memclass) const type memclass

/* VAR macro for the declaration and definition of variables
   vartype        type of the variable
   memclass    classification of the variable itself
 */
# define VAR(vartype, memclass) vartype memclass


/*To be used for code*/
#define STRS_CODE

/* To be used for all global or static variables that have at least one of the
following properties:
-accessed bitwise
-frequently used
-high number of accesses in source code*/
#define STRS_VAR_FAST

/*To be used for global or static variables that are initialized after every reset.*/
#define STRS_VAR

/*To be used for global or static constants*/
#define STRS_CONST


#endif /* COMPILER_CFG_H_ */
