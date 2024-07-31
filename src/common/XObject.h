#ifndef XOBJECT_H
#define XOBJECT_H

#include <iostream>

class XObject {
    public:
        bool m_verbose = false;
        bool m_valid = false;
        bool m_reuse = false;
		const char * m_name;  // or static?

        void set_verbose(bool verbose){
            this->m_verbose = verbose;
        }

		bool get_verbose(){
            return this->m_verbose();
        }

        void set_name(const char * name){
            this->m_name = name;
        }

		bool get_verbose(){
            return this->m_verbose();
        }

        virtual void init() = 0; // constructor
        virtual void release() = 0; // deconstruct

        virtual void print_settings() = 0;
}
#endif //XOBJECT_H
