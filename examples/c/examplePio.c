#include <stdio.h>
#include <stdlib.h>

typedef struct examplePioClass
{
    int someThing;
    struct examplePioClass* (*init)  (struct examplePioClass*);
    struct examplePioClass* (*write) (struct examplePioClass*);
    struct examplePioClass* (*read)  (struct examplePioClass*);
    
} examplePioClass;

struct examplePioClass* epc_init( struct examplePioClass* this )
{
    /* implementation of init */
    printf(" examplePioClass::init %d\n",this->someThing);
    return this;
}

struct examplePioClass* epc_write( struct examplePioClass* this )
{
    /* implementation of write */
    printf(" examplePioClass::write %d\n",this->someThing);
    return this;
}

struct examplePioClass* epc_read( struct examplePioClass* this )
{
    /* implementation of read */
    printf(" examplePioClass::read %d\n",this->someThing);
    return this;
}

struct examplePioClass* epc_new()
{
    /* assign function pointers to impls */
    printf(" examplePioClass::new \n");
    struct examplePioClass* this = malloc((sizeof(struct examplePioClass)));
    this->init = epc_init;
    this->write = epc_write;
    this->read = epc_read;
    this->someThing = 10;
    return this;
}

int main(int argc, const char* argv[])
{
    struct examplePioClass* pioExInst = epc_new();
    
    pioExInst->init(pioExInst);
    pioExInst->write(pioExInst);
    pioExInst->read(pioExInst);

    free(pioExInst);
    return 0;

}