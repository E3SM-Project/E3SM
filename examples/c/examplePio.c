#include <stdio.h>
#include <stdlib.h>

typedef struct examplePioClass
{
    void (*init)  (struct examplePioClass*);
    void (*write) (struct examplePioClass*);
    void (*read)  (struct examplePioClass*);
    
} examplePioClass;

void epc_init( struct examplePioClass* this )
{
    /* implementation of init */
    printf(" examplePioClass::init \n");
    //return this;
}

void epc_write( struct examplePioClass* this )
{
    /* implementation of write */
    printf(" examplePioClass::write \n");
    //return this;
}

void epc_read( struct examplePioClass* this )
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
    
    pioExInst->init(pioExInst);
    pioExInst->write(pioExInst);
    pioExInst->read(pioExInst);

    free(pioExInst);
    return 0;

}