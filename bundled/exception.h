/**
 * @file exception.h
 * @brief  一些简单的诊断函数
 * @author Yang Fanyi <yfy__92@126.com>
 * @version 
 */


#ifndef __EXCEPTION_H_
#define __EXCEPTION_H_

#include <cassert>
#include <iostream>

#ifdef DEBUG
#define Assert(cond, error) \
{ \
  if(!(cond)){\
    std::cerr << error << std::endl; \
    assert(cond); \
  }}
#else
#define Assert(cond, error) {}
#endif

#define AssertExc(cond, error) \
{ \
  if(!(cond)){\
    std::cerr << "An error occurred in line <" << __LINE__ \
    << "> of file <" << __FILE__ << ">" << std::endl << \
    "The violated condition is : " << #cond << std::endl << \
    "The description is : " << error << std::endl; \
    std::abort();\
  }}

#define AssertCond(cond) \
{ \
  if(!(cond)){\
    std::cerr << "An error occurred in line <" << __LINE__ \
    << "> of file <" << __FILE__ << ">" << std::endl << \
    "The violated condition is : " << #cond << std::endl << \
    "The description is : " << #cond << std::endl; \
    std::abort();\
  }}

#define AssertExp0(cond) AssertCond(cond)

#define AssertExp1(cond, arg1) \
{ \
    if(!(cond)){\
    std::cerr << "An error occurred in line <" << __LINE__ \
    << "> of file <" << __FILE__ << ">" << std::endl << \
    "The violated condition is : " << #cond << std::endl << \
    "The description is : " << #cond << std::endl << \
    "The arg1 is : " << arg1 << std::endl; \
    std::abort();\
  }}

#define AssertExp2(cond, arg1, arg2) \
{ \
    if(!(cond)){\
    std::cerr << "An error occurred in line <" << __LINE__ \
    << "> of file <" << __FILE__ << ">" << std::endl << \
    "The violated condition is : " << #cond << std::endl << \
    "The description is : " << #cond << std::endl << \
    "The arg1 is : " << arg1 << std::endl << \
    "The arg2 is : " << arg2 << std::endl; \
    std::abort();\
  }}

#define AssertExp3(cond, arg1, arg2, arg3) \
{ \
    if(!(cond)){\
    std::cerr << "An error occurred in line <" << __LINE__ \
    << "> of file <" << __FILE__ << ">" << std::endl << \
    "The violated condition is : " << #cond << std::endl << \
    "The description is : " << #cond << std::endl << \
    "The arg1 is : " << arg1 << std::endl << \
    "The arg2 is : " << arg2 << std::endl << \
    "The arg3 is : " << arg3 << std::endl; \
    std::abort();\
  }}

#endif
