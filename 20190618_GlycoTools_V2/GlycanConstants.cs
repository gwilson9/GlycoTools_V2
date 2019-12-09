﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _20190618_GlycoTools_V2
{
    public class GlycanConstants
    {
        public readonly List<string> glycans = "HexNAc(1);HexNAc(1)Fuc(1);HexNAc(2);HexNAc(2)Fuc(1);HexNAc(2)Hex(1);HexNAc(2)Hex(1)Fuc(1);HexNAc(2)Hex(2);HexNAc(2)Hex(2)Fuc(1);HexNAc(2)Hex(3);HexNAc(2)Hex(3)Fuc(1);HexNAc(2)Hex(4);HexNAc(2)Hex(4)Fuc(1);HexNAc(2)Hex(5);HexNAc(2)Hex(5)Fuc(1);HexNAc(2)Hex(6);HexNAc(2)Hex(6)Fuc(1);HexNAc(2)Hex(6)Phospho(1);HexNAc(2)Hex(7);HexNAc(2)Hex(7)Fuc(1);HexNAc(2)Hex(8);HexNAc(2)Hex(9);HexNAc(2)Hex(10);HexNAc(2)Hex(11);HexNAc(3)Hex(3);HexNAc(3)Hex(3)Fuc(1);HexNAc(3)Hex(4);HexNAc(3)Hex(4)Fuc(1);HexNAc(3)Hex(4)NeuAc(1);HexNAc(3)Hex(4)Fuc(1)NeuAc(1);HexNAc(3)Hex(4)Fuc(2)NeuAc(1);HexNAc(3)Hex(5);HexNAc(3)Hex(5)Fuc(1);HexNAc(3)Hex(5)NeuAc(1);HexNAc(3)Hex(5)Fuc(1)NeuAc(1);HexNAc(3)Hex(5)Fuc(1)NeuAc(2);HexNAc(3)Hex(6);HexNAc(3)Hex(6)Fuc(1)NeuAc(1);HexNAc(3)Hex(6)NeuAc(1);HexNAc(4)Hex(3);HexNAc(4)Hex(3)Fuc(1);HexNAc(4)Hex(4);HexNAc(4)Hex(4)Fuc(1);HexNAc(4)Hex(4)NeuAc(1);HexNAc(4)Hex(5);HexNAc(4)Hex(5)Fuc(1);HexNAc(4)Hex(6)Fuc(1);HexNAc(4)Hex(6)Fuc(2);HexNAc(4)Hex(5)NeuAc(1);HexNAc(4)Hex(5)NeuAc(2);HexNAc(4)Hex(5)Fuc(1)NeuAc(1);HexNAc(4)Hex(5)Fuc(2)NeuAc(1);HexNAc(4)Hex(5)Fuc(3)NeuAc(1);HexNAc(4)Hex(5)Fuc(1)NeuAc(2);HexNAc(4)Hex(5)Fuc(2)NeuAc(2);HexNAc(4)Hex(5)Fuc(3)NeuAc(2);HexNAc(5)Hex(3);HexNAc(5)Hex(3)Fuc(1);HexNAc(5)Hex(3)Fuc(2);HexNAc(5)Hex(4);HexNAc(5)Hex(4)Fuc(1);HexNAc(5)Hex(4)Fuc(2);HexNAc(5)Hex(5);HexNAc(5)Hex(5)Fuc(1);HexNAc(5)Hex(5)Fuc(2);HexNAc(5)Hex(5)Fuc(3);HexNAc(5)Hex(5)NeuAc(1);HexNAc(5)Hex(5)NeuAc(2);HexNAc(5)Hex(5)Fuc(1)NeuAc(1);HexNAc(5)Hex(5)Fuc(2)NeuAc(1);HexNAc(5)Hex(5)Fuc(1)NeuAc(2);HexNAc(5)Hex(5)Fuc(2)NeuAc(2);HexNAc(5)Hex(6);HexNAc(5)Hex(6)Fuc(1);HexNAc(5)Hex(6)Fuc(2);HexNAc(5)Hex(6)Fuc(3);HexNAc(5)Hex(6)Fuc(4);HexNAc(5)Hex(6)NeuAc(1);HexNAc(5)Hex(6)Fuc(1)NeuAc(1);HexNAc(5)Hex(6)Fuc(2)NeuAc(1);HexNAc(5)Hex(6)Fuc(3)NeuAc(1);HexNAc(5)Hex(6)Fuc(4)NeuAc(1);HexNAc(5)Hex(6)NeuAc(2);HexNAc(5)Hex(6)Fuc(1)NeuAc(2);HexNAc(5)Hex(6)Fuc(2)NeuAc(2);HexNAc(5)Hex(6)Fuc(3)NeuAc(2);HexNAc(5)Hex(6)Fuc(4)NeuAc(2);HexNAc(5)Hex(6)NeuAc(3);HexNAc(5)Hex(6)Fuc(1)NeuAc(3);HexNAc(5)Hex(6)Fuc(2)NeuAc(3);HexNAc(5)Hex(6)Fuc(3)NeuAc(3);HexNAc(5)Hex(6)Fuc(4)NeuAc(3);HexNAc(6)Hex(3);HexNAc(6)Hex(3)Fuc(1);HexNAc(6)Hex(3)Fuc(2);HexNAc(6)Hex(4);HexNAc(6)Hex(4)Fuc(1);HexNAc(6)Hex(4)Fuc(2);HexNAc(6)Hex(5);HexNAc(6)Hex(5)Fuc(1);HexNAc(6)Hex(5)Fuc(2);HexNAc(6)Hex(5)Fuc(2)NeuAc(1);HexNAc(6)Hex(5)Fuc(3);HexNAc(6)Hex(6);HexNAc(6)Hex(6)Fuc(1);HexNAc(6)Hex(6)Fuc(2);HexNAc(6)Hex(6)NeuAc(1);HexNAc(6)Hex(6)Fuc(1)NeuAc(1);HexNAc(6)Hex(6)Fuc(2)NeuAc(1);HexNAc(6)Hex(6)NeuAc(2);HexNAc(6)Hex(6)Fuc(1)NeuAc(2);HexNAc(6)Hex(6)Fuc(2)NeuAc(2);HexNAc(6)Hex(6)NeuAc(3);HexNAc(6)Hex(6)Fuc(1)NeuAc(3);HexNAc(6)Hex(6)Fuc(2)NeuAc(3);HexNAc(6)Hex(6)Fuc(3)NeuAc(3);HexNAc(6)Hex(7)Fuc(1);HexNAc(6)Hex(7)Fuc(2);HexNAc(6)Hex(7)Fuc(3);HexNAc(6)Hex(7)NeuAc(1);HexNAc(6)Hex(7)Fuc(1)NeuAc(1);HexNAc(6)Hex(7)Fuc(2)NeuAc(1);HexNAc(6)Hex(7)Fuc(3)NeuAc(1);HexNAc(6)Hex(7)NeuAc(2);HexNAc(6)Hex(7)Fuc(1)NeuAc(2);HexNAc(6)Hex(7)Fuc(2)NeuAc(2);HexNAc(6)Hex(7)Fuc(3)NeuAc(2);HexNAc(6)Hex(7)NeuAc(3);HexNAc(6)Hex(7)Fuc(1)NeuAc(3);HexNAc(6)Hex(7)Fuc(2)NeuAc(3);HexNAc(6)Hex(7)Fuc(3)NeuAc(3);HexNAc(6)Hex(7)NeuAc(4);HexNAc(6)Hex(7)Fuc(1)NeuAc(4);HexNAc(6)Hex(7)Fuc(2)NeuAc(4);HexNAc(6)Hex(7)Fuc(3)NeuAc(4);HexNAc(4)Hex(3)NeuAc(1);HexNAc(3)Hex(6)Fuc(1);HexNAc(4)Hex(6);HexNAc(4)Hex(4)Fuc(1)NeuAc(1);HexNAc(7)Hex(3);HexNAc(5)Hex(3)Fuc(1)NeuAc(1);HexNAc(4)Hex(7);HexNAc(5)Hex(4)NeuAc(1);HexNAc(7)Hex(3)Fuc(1);HexNAc(7)Hex(4);HexNAc(4)Hex(6)NeuAc(1);HexNAc(4)Hex(7)Fuc(1);HexNAc(5)Hex(4)Fuc(1)NeuAc(1);HexNAc(8)Hex(3);HexNAc(6)Hex(3)Fuc(1)NeuAc(1);HexNAc(5)Hex(7);HexNAc(6)Hex(4)NeuAc(1);HexNAc(7)Hex(4)Fuc(1);HexNAc(4)Hex(6)Fuc(1)NeuAc(1);HexNAc(4)Hex(7)NeuAc(1);HexNAc(5)Hex(4)NeuAc(2);HexNAc(8)Hex(3)Fuc(1);HexNAc(8)Hex(4);HexNAc(5)Hex(7)Fuc(1);HexNAc(5)Hex(8);HexNAc(9)Hex(3);HexNAc(2)Hex(12);HexNAc(6)Hex(7);HexNAc(5)Hex(4)Fuc(1)NeuAc(2);HexNAc(7)Hex(6);HexNAc(6)Hex(3)Fuc(1)NeuAc(2);HexNAc(8)Hex(5);HexNAc(5)Hex(8)Fuc(1);HexNAc(9)Hex(3)Fuc(1);HexNAc(9)Hex(4);HexNAc(7)Hex(6)Fuc(1);HexNAc(7)Hex(7);HexNAc(8)Hex(5)Fuc(1);HexNAc(5)Hex(7)Fuc(1)NeuAc(1);HexNAc(8)Hex(6);HexNAc(5)Hex(9)Fuc(1);HexNAc(9)Hex(4)Fuc(1);HexNAc(6)Hex(9);HexNAc(7)Hex(7)Fuc(1);HexNAc(7)Hex(8);HexNAc(8)Hex(7);HexNAc(9)Hex(6);HexNAc(6)Hex(8)NeuAc(1);HexNAc(7)Hex(8)Fuc(1);HexNAc(5)Hex(7)Fuc(1)NeuAc(2);HexNAc(8)Hex(8);HexNAc(9)Hex(6)Fuc(1);HexNAc(6)Hex(8)Fuc(1)NeuAc(1);HexNAc(7)Hex(8)NeuAc(1);HexNAc(6)Hex(5)Fuc(1)NeuAc(2);HexNAc(6)Hex(5)Fuc(1)NeuAc(1);HexNAc(6)Hex(5)Fuc(1)NeuAc(3);HexNAc(8)Hex(8)Fuc(1);HexNAc(8)Hex(9);HexNAc(6)Hex(11)Fuc(1);HexNAc(10)Hex(7);HexNAc(8)Hex(9)Fuc(1);HexNAc(6)Hex(10)Fuc(1)NeuAc(1);HexNAc(7)Hex(7)Fuc(1)NeuAc(2);HexNAc(6)Hex(9)Fuc(1)NeuAc(2);HexNAc(9)Hex(9)Fuc(1);HexNAc(7)Hex(8)Fuc(1)NeuAc(1);HexNAc(7)Hex(8)Fuc(1)NeuAc(2);HexNAc(9)Hex(10);HexNAc(7)Hex(7)Fuc(1)NeuAc(3);HexNAc(9)Hex(10)Fuc(1);HexNAc(7)Hex(6)Fuc(1)NeuAc(4);HexNAc(7)Hex(8)Fuc(1)NeuAc(3);HexNAc(10)Hex(10)Fuc(1);HexNAc(7)Hex(7)Fuc(1)NeuAc(4);HexNAc(7)Hex(8)Fuc(1)NeuAc(4);HexNAc(8)Hex(9)Fuc(1)NeuAc(3);HexNAc(11)Hex(11)NeuAc(1);HexNAc(8)Hex(9)Fuc(1)NeuAc(4);HexNAc(9)Hex(10)Fuc(1)NeuAc(4);HexNAc(2)Hex(13);HexNAc(2)Hex(14);HexNAc(2)Hex(15);HexNAc(2)Hex(16);HexNAc(2)Hex(17);HexNAc(2)Hex(18);HexNAc(1)Fuc(2);HexNAc(2)Fuc(2);HexNAc(2)Hex(1)Fuc(2);HexNAc(2)Hex(2)Fuc(2);HexNAc(2)Hex(3)Fuc(2);HexNAc(2)Hex(4)Fuc(2);HexNAc(2)Hex(5)Fuc(2);HexNAc(3)Hex(3)Fuc(2);HexNAc(3)Hex(4)Fuc(1);HexNAc(3)Hex(4);HexNAc(3)Hex(4)Fuc(2);HexNAc(3)Hex(4)NeuGc(1);HexNAc(4)Hex(3)Fuc(2);HexNAc(4)Hex(3)NeuGc(1);HexNAc(3)Hex(4)Fuc(1)NeuGc(1);HexNAc(3)Hex(5)NeuGc(1);HexNAc(4)Hex(3)Fuc(3);HexNAc(4)Hex(4)Fuc(2);HexNAc(4)Hex(4)NeuGc(1);HexNAc(3)Hex(6)NeuGc(1);HexNAc(4)Hex(4)Fuc(1)NeuGc(1);HexNAc(4)Hex(5)Fuc(2);HexNAc(4)Hex(5)NeuGc(1);HexNAc(5)Hex(3)Fuc(1)NeuGc(1);HexNAc(5)Hex(4)NeuGc(1);HexNAc(3)Hex(6)Fuc(1)NeuGc(1);HexNAc(4)Hex(4)Fuc(2)NeuAc(1);HexNAc(4)Hex(4)Fuc(2)NeuGc(1);HexNAc(4)Hex(5)Fuc(3);HexNAc(4)Hex(5)Fuc(1)NeuGc(1);HexNAc(4)Hex(6)NeuGc(1);HexNAc(5)Hex(4)Fuc(1)NeuGc(1);HexNAc(5)Hex(5)NeuGc(1);HexNAc(6)Hex(3)Fuc(3);HexNAc(6)Hex(3)Fuc(1)NeuGc(1);HexNAc(4)Hex(5)Fuc(4);HexNAc(4)Hex(5)NeuAc(1)NeuGc(1);HexNAc(4)Hex(5)Fuc(2)NeuGc(1);HexNAc(4)Hex(6)Fuc(3);HexNAc(4)Hex(5)NeuGc(2);HexNAc(4)Hex(6)Fuc(1)NeuGc(1);HexNAc(4)Hex(7)Fuc(2);HexNAc(5)Hex(4)Fuc(2)NeuAc(1);HexNAc(5)Hex(4)NeuAc(1)NeuGc(1);HexNAc(5)Hex(4)Fuc(2)NeuGc(1);HexNAc(5)Hex(4)NeuGc(2);HexNAc(5)Hex(5)Fuc(1)NeuGc(1);HexNAc(6)Hex(3)Fuc(2)NeuAc(1);HexNAc(5)Hex(6)NeuGc(1);HexNAc(6)Hex(3)Fuc(2)NeuGc(1);HexNAc(7)Hex(4)Fuc(2);HexNAc(4)Hex(5)Fuc(1)NeuAc(1)NeuGc(1);HexNAc(4)Hex(5)Fuc(1)NeuGc(2);HexNAc(5)Hex(4)Fuc(1)NeuAc(1)NeuGc(1);HexNAc(5)Hex(4)Fuc(1)NeuGc(2);HexNAc(5)Hex(5)NeuAc(1)NeuGc(1);HexNAc(5)Hex(5)NeuGc(2);HexNAc(5)Hex(6)Fuc(1)NeuGc(1);HexNAc(6)Hex(3)Fuc(1)NeuAc(1)NeuGc(1);HexNAc(6)Hex(3)Fuc(1)NeuGc(2);HexNAc(5)Hex(5)Fuc(1)NeuAc(1)NeuGc(1);HexNAc(5)Hex(5)Fuc(1)NeuGc(2);HexNAc(5)Hex(6)NeuAc(1)NeuGc(1);HexNAc(5)Hex(6)Fuc(2)NeuGc(1);HexNAc(5)Hex(6)NeuGc(2);HexNAc(6)Hex(6)Fuc(3);HexNAc(4)Hex(5)Fuc(3)NeuAc(1)NeuGc(1);HexNAc(6)Hex(7)NeuGc(1);HexNAc(4)Hex(5)Fuc(3)NeuGc(2);HexNAc(5)Hex(6)Fuc(1)NeuAc(1)NeuGc(1);HexNAc(5)Hex(6)Fuc(3)NeuGc(1);HexNAc(5)Hex(6)Fuc(1)NeuGc(2);HexNAc(6)Hex(6)NeuAc(1)NeuGc(1);HexNAc(6)Hex(6)NeuGc(2);HexNAc(6)Hex(7)Fuc(1)NeuGc(1);HexNAc(5)Hex(6)NeuAc(2)NeuGc(1);HexNAc(5)Hex(6)NeuAc(1)NeuGc(2);HexNAc(5)Hex(7)Fuc(1)NeuAc(1)NeuGc(1);HexNAc(5)Hex(8)Fuc(4);HexNAc(5)Hex(6)NeuGc(3);HexNAc(5)Hex(7)Fuc(1)NeuGc(2);HexNAc(6)Hex(6)Fuc(1)NeuAc(1)NeuGc(1);HexNAc(6)Hex(7)Fuc(4);HexNAc(6)Hex(7)NeuAc(1)NeuGc(1);HexNAc(6)Hex(6)Fuc(1)NeuGc(2);HexNAc(6)Hex(7)NeuGc(2);HexNAc(5)Hex(6)Fuc(1)NeuAc(2)NeuGc(1);HexNAc(7)Hex(8)NeuGc(1);HexNAc(5)Hex(6)Fuc(1)NeuAc(1)NeuGc(2);HexNAc(5)Hex(6)Fuc(1)NeuGc(3);HexNAc(6)Hex(5)Fuc(1)NeuAc(1)NeuGc(1);HexNAc(6)Hex(5)Fuc(1)NeuAc(2)NeuGc(1);HexNAc(6)Hex(5)Fuc(1)NeuAc(1)NeuGc(2);HexNAc(6)Hex(7)Fuc(5);HexNAc(6)Hex(5)Fuc(1)NeuGc(3);HexNAc(6)Hex(6)Fuc(1)NeuAc(2)NeuGc(1);HexNAc(6)Hex(7)Fuc(4)NeuAc(1);HexNAc(6)Hex(7)NeuAc(2)NeuGc(1);HexNAc(6)Hex(6)Fuc(1)NeuAc(1)NeuGc(2);HexNAc(6)Hex(7)Fuc(4)NeuGc(1);HexNAc(6)Hex(6)Fuc(1)NeuGc(3);HexNAc(6)Hex(10)Fuc(1)NeuGc(1);HexNAc(6)Hex(9)Fuc(1)NeuAc(1)NeuGc(1);HexNAc(7)Hex(8)Fuc(1)NeuAc(1)NeuGc(1);HexNAc(7)Hex(8)Fuc(1)NeuGc(2);HexNAc(6)Hex(7)NeuAc(3)NeuGc(1);HexNAc(6)Hex(7)NeuAc(2)NeuGc(2);HexNAc(7)Hex(7)Fuc(1)NeuAc(2)NeuGc(1);HexNAc(7)Hex(7)Fuc(1)NeuAc(1)NeuGc(2);HexNAc(6)Hex(7)Fuc(1)NeuAc(3)NeuGc(1);HexNAc(7)Hex(6)Fuc(1)NeuAc(3)NeuGc(1);HexNAc(7)Hex(8)Fuc(1)NeuAc(2)NeuGc(1);HexNAc(7)Hex(8)Fuc(1)NeuAc(1)NeuGc(2);HexNAc(7)Hex(8)Fuc(1)NeuGc(3);HexNAc(8)Hex(9)Fuc(1)NeuAc(2)NeuGc(1);HexNAc(8)Hex(9)Fuc(1)NeuAc(1)NeuGc(2);HexNAc(11)Hex(11)NeuGc(1);HexNAc(8)Hex(9)Fuc(1)NeuAc(2);HexNAc(9)Hex(10)Fuc(1)NeuAc(3);HexNAc(4)Hex(4)Fuc(1)Pent(1);HexNAc(4)Hex(5)NeuAc(1);HexNAc(4)Hex(5)NeuAc(2);HexNAc(4)Hex(5)NeuAc(1)NeuGc(1);HexNAc(4)Hex(5)NeuGc(2);HexNAc(4)Hex(5)Fuc(1)NeuGc(2);HexNAc(3)Hex(2)Fuc(1);HexNAc(3)Hex(2);HexNAc(2)Hex(6)Fuc(2);HexNAc(3)Hex(5)Fuc(2);HexNAc(2)Hex(3)Pent(1);HexNAc(3)Hex(3)Pent(1);HexNAc(4)Hex(3)Pent(1);HexNAc(2)Hex(3)Fuc(1)Pent(1);HexNAc(3)Hex(3)Fuc(1)Pent(1);HexNAc(4)Hex(3)Fuc(1)Pent(1);HexNAc(3)Hex(3)Fuc(2)Pent(1);HexNAc(3)Hex(4)Fuc(1)Pent(1);HexNAc(3)Hex(4)Fuc(2)Pent(1);HexNAc(4)Hex(3)Fuc(2)Pent(1);HexNAc(4)Hex(4)Fuc(2)Pent(1);HexNAc(4)Hex(5)Fuc(1)Pent(1);HexNAc(4)Hex(5)Fuc(2)Pent(1);HexNAc(4)Hex(4)Fuc(3)Pent(1);HexNAc(4)Hex(3)Fuc(3)Pent(1);HexNAc(4)Hex(5)Fuc(3)Pent(1);HexNAc(3)Hex(9);HexNAc(2)Hex(4)Fuc(1)Pent(1);HexNAc(6)Hex(4)NeuAc(2);HexNAc(4)Hex(5)NeuAc(3);HexNAc(4)Hex(5)NeuAc(2)NeuGc(1);HexNAc(4)Hex(5)NeuAc(1)NeuGc(2);HexNAc(4)Hex(5)NeuAc(4);HexNAc(5)Hex(6)NeuAc(4);HexNAc(5)Hex(6)Fuc(1)NeuAc(4);HexNAc(5)Hex(6)Fuc(1)NeuAc(5);HexNAc(6)Hex(7)NeuAc(5);HexNAc(6)Hex(7)Fuc(1)NeuAc(5);Hex(1);Hex(2);Hex(2)Phospho(1);Hex(3);Hex(4);Hex(5);Hex(6);Hex(7);Hex(8);Hex(9);Hex(10);Hex(11);Hex(12);HexNAc(1)Hex(1);HexNAc(3)Hex(1);HexNAc(1)Hex(1)NeuAc(1);HexNAc(1)Hex(1)NeuAc(2);HexNAc(2)Hex(2)NeuAc(1);HexNAc(3)Hex(1)NeuAc(2);HexNAc(2)Hex(2)Fuc(1)NeuAc(1);HexNAc(2)Hex(2)NeuAc(2);HexNAc(1)Hex(1)NeuAc(3);Pent(1);Pent(2);Pent(3);Pent(4);Pent(5);HexNAc(1)Hex(1)Fuc(1);HexNAc(3);HexNAc(1)Hex(2)Fuc(1);HexNAc(3)Fuc(1);HexNAc(1)Hex(1)Fuc(1)NeuAc(1);HexNAc(2)Hex(1)NeuAc(1);HexNAc(3)Hex(1)Fuc(1);HexNAc(4)Hex(1);HexNAc(2)Hex(1)Fuc(1)NeuAc(1);HexNAc(3)Hex(1)NeuAc(1);HexNAc(3)Hex(1)Fuc(2);HexNAc(4)Hex(1)Fuc(1);HexNAc(4)Hex(2);HexNAc(2)Hex(1)NeuAc(2);HexNAc(2)Hex(2)Fuc(3);HexNAc(2)Hex(3)NeuAc(1);HexNAc(3)Hex(1)Fuc(1)NeuAc(1);HexNAc(3)Hex(2)NeuAc(1);HexNAc(3)Hex(2)Fuc(2);HexNAc(2)Hex(2)Fuc(2)NeuAc(1);HexNAc(3)Hex(2)Fuc(1)NeuAc(1);HexNAc(3)Hex(2)Fuc(3);HexNAc(3)Hex(3)NeuAc(1);HexNAc(4)Hex(2)Fuc(2);HexNAc(2)Hex(2)Fuc(1)NeuAc(2);HexNAc(3)Hex(2)NeuAc(2);HexNAc(3)Hex(2)Fuc(2)NeuAc(1);HexNAc(3)Hex(2)Fuc(4);HexNAc(3)Hex(3)Fuc(1)NeuAc(1);HexNAc(3)Hex(3)Fuc(3);HexNAc(3)Hex(3)NeuAc(2);HexNAc(3)Hex(3)Fuc(1)NeuAc(2);HexNAc(4)Hex(4)Fuc(3);HexNAc(3)Hex(3)Fuc(2)NeuAc(2);HexNAc(4)Hex(4)Fuc(3)NeuAc(1);HexNAc(1)Hex(1)NeuGc(1);HexNAc(1)Hex(1)Fuc(1)NeuGc(1);HexNAc(2)Hex(1)NeuGc(1);HexNAc(1)Hex(1)NeuGc(1)NeuAc(1);HexNAc(1)Hex(1)NeuGc(2);HexNAc(2)Hex(1)Fuc(1)NeuGc(1);HexNAc(2)Hex(1)Fuc(2)NeuGc(1);HexNAc(3)Hex(1)Fuc(1)NeuGc(1);Fuc(1)".Split(';').ToList();

    }
}
