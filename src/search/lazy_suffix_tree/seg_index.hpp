#pragma once

namespace idx {
    struct segmented_long {
        uint8_t segment;
        uint8_t index_0;
        uint8_t index_1;
        uint8_t index_2;
        uint8_t index_3;
    };

    segmented_long deconstruct_index(size_t real_index) {
      segmented_long return_struct;
      uint8_t segment_ = real_index >> 32;
      unsigned int index_ = real_index & 0x0FFFFFFFF;
      int high16 = index_ >> 16;
      int low16 = index_ & 0xFFFF;
      uint8_t index_0 = low16 & 0xFF;
      uint8_t index_1 = low16 >> 8;
      uint8_t index_2 = high16 & 0xFF;
      uint8_t index_3 = high16 >> 8;
      return_struct.index_0 = index_0;
      return_struct.index_1 = index_1;
      return_struct.index_2 = index_2;
      return_struct.index_3 = index_3;
      return_struct.segment = segment_;
      return return_struct;
    }

    size_t reconstruct_index(segmented_long seg_i) {
      size_t base = seg_i.segment;
      base = base << 32;
      int high16 = (seg_i.index_3 << 8) | seg_i.index_2;
      int low16 = (seg_i.index_1 << 8) | seg_i.index_0;
      int index = (high16 << 16) | (low16 & 0xFFFF);
      return (base + index);
    }
}
