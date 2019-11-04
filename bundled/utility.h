/**
 * @file utility.h
 * @brief 一些辅助工具
 * @author Yang Fanyi <yfy__92@126.com>
 * @version 
 */


#ifndef __UTILITY_H_
#define __UTILITY_H_

#include "exception.h"

namespace Utilities{
  template <typename Iterator, typename T, typename Comp>
  inline Iterator
  lower_bound(Iterator first, Iterator last, 
      const T& val, const Comp comp){
    Assert(last - first >= 0, "the iterators do not exist");

    unsigned int len = static_cast<unsigned int>(last - first);
    if(len == 0) return first;

    while(true){
      if(len < 8){
        switch(len){
          case 7:
            if(!comp(*first, val))
              return first;
            ++first;
          case 6:
            if(!comp(*first, val))
              return first;
            ++first;
          case 5:
            if(!comp(*first, val))
              return first;
            ++first;
          case 4:
            if(!comp(*first, val))
              return first;
            ++first;
          case 3:
            if(!comp(*first, val))
              return first;
            ++first;
          case 2:
            if(!comp(*first, val))
              return first;
            ++first;
          case 1:
            if(!comp(*first, val))
              return first;
            return first + 1;
          default:
            std::abort();
        }
      }

      const unsigned int half = len >> 1;
      const Iterator middle = first + half;
      if(comp(*middle, val)){
        first = middle + 1;
        len -= half + 1;
      }
      else
        len = half;
    }
  }

  template <typename Iterator, typename T>
  inline Iterator
  lower_bound(Iterator first, Iterator last, const T& val){
    return Utilities::lower_bound(first, last, val, std::less<T>());
  }

};

#endif
