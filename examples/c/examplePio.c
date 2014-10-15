#include <stdio.h>
#include <stdlib.h>

typedef struct examplePioClass
{
    void (*init)  ();
    void (*write) ();
    void (*read)  ();
    
} examplePioClass;

void epc_init( )
{
    /* implementation of init */
    printf(" examplePioClass::init \n");
    //return this;
}

void epc_write( )
{
    /* implementation of write */
    printf(" examplePioClass::write \n");
    //return this;
}

void epc_read( )
{
    /* implementation of read */
    printf(" examplePioClass::read \n");
    //return this;
}

struct examplePioClass* epc_new()
{
    /* assign function pointers to impls */
    printf(" examplePioClass::new \n");
    struct examplePioClass* this = malloc((sizeof(struct examplePioClass)));
    this->init = epc_init;
    this->write = epc_write;
    this->read = epc_read;
    return this;
}

int main(int argc, const char* argv[])
{
    struct examplePioClass* pioExInst = epc_new();
    
    pioExInst->init();
    pioExInst->write();
    pioExInst->read();

    free(pioExInst);
    return 0;

}