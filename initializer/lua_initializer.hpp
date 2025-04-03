#ifndef UGPLUGIN_XBRAIDFORUG4_INITIALIZER_LUA_INITIALIZER_HPP
#define UGPLUGIN_XBRAIDFORUG4_INITIALIZER_LUA_INITIALIZER_HPP

// todo implement class


void setGeneratorComponent(const char *cmp) {
    this->cmp_ = cmp;
}
#ifdef UG_FOR_LUA
    void setVectorGenerator(const char *fctName) {
        setVectorGenerator(ug::LuaUserDataFactory<double, TDomain::dim>::create(fctName));
    }

void setVectorGenerator(ug::LuaFunctionHandle fct) {
        setVectorGenerator(make_sp(new ug::LuaUserData<double, TDomain::dim>(fct)));
    }
#endif




#endif
