#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 22:19:06 2026

@author: eric
"""

from enum import IntFlag
    
class CellType(IntFlag):
    # flag during process (be careful in binary test)
    TOBEPROCESSED = 1
    
    nanshift = 2
    # Out of domain points, distinguish points in holes and outside
    INNaN = nanshift # df.NOTIN for points that are in holes
    OUT = nanshift << 1 # for points that are outside the outer boundary
    NaN = INNaN | OUT # No value points

    # Point/line flag
    POINT = OUT << 1
    LINE = POINT << 1
    POLYGON = LINE << 1
    
    pointshift = POLYGON << 1
    # Regular point
    REGULAR = POINT | pointshift 
    
    # Critical points
    SADDLE_FLAG = pointshift << 1
    PIT_FLAG = pointshift << 2
    PEAK_FLAG = pointshift << 3
    TJUNCTION_FLAG = pointshift << 4
    RJUNCTION_FLAG = pointshift << 5
    CULVERTNODE_FLAG = pointshift << 6
    SADDLE = POINT | SADDLE_FLAG
    PIT = POINT | PIT_FLAG
    PEAK = POINT | PEAK_FLAG
    TJUNCTION = POINT | TJUNCTION_FLAG
    RJUNCTION = POINT | RJUNCTION_FLAG
    CULVERTNODE = POINT | CULVERTNODE_FLAG
    COURSENODE_FLAG = SADDLE_FLAG | PIT_FLAG | TJUNCTION_FLAG | CULVERTNODE_FLAG
    DRYNODE_FLAG = PEAK_FLAG | RJUNCTION_FLAG
    
#    fluence = pointshift << 6 # à voir si nécessaire de les garder
    # Fluences
#    TRANSFLUENCE = POINT | fluence
#    SINK = POINT | fluence << 1
#    SPRING = POINT | fluence << 2
#    CONFLUENCE = POINT | fluence << 3
#    BIFURCATION = POINT | fluence << 4    
    
    lineshift = pointshift
    # Main types of line
    DIVIDE_FLAG = lineshift
    DIVIDE = LINE | DIVIDE_FLAG
    WATERCOURSE_FLAG = lineshift << 1
    WATERCOURSE = LINE | WATERCOURSE_FLAG

    # Flags under divide
    RIDGE_FLAG = lineshift << 2
    # Flags under watercourse
    SURFACE_FLAG = lineshift << 3
    UNDERGROUND_FLAG = lineshift << 4
    # Flags at third level
    THALWEG_FLAG = lineshift << 5
    CULVERT_FLAG = lineshift << 6

    # Subtypes built with the flags
    RIDGE = DIVIDE | RIDGE_FLAG
    SURFACE = WATERCOURSE | SURFACE_FLAG
    UNDERGROUND = WATERCOURSE | UNDERGROUND_FLAG
    THALWEG = SURFACE     | THALWEG_FLAG
    CULVERT = UNDERGROUND | CULVERT_FLAG
    
    @staticmethod
    def isSaddle(flagvalue: "CellType") -> bool:
        return flagvalue & CellType.SADDLE_FLAG and flagvalue & CellType.POINT 

    @staticmethod
    def isThalweg(flagvalue: "CellType") -> bool:
        return flagvalue & CellType.THALWEG_FLAG and flagvalue & CellType.LINE

    @staticmethod
    def isSurfaceWatercourse(flagvalue: "CellType") -> bool:
        return flagvalue & CellType.SURFACE_FLAG and flagvalue & CellType.LINE

    @staticmethod
    def isTJunction(flagvalue: "CellType") -> bool:
        return flagvalue & CellType.TJUNCTION_FLAG and flagvalue & CellType.POINT

    @staticmethod
    def isRidge(flagvalue: "CellType") -> bool:
        return flagvalue & CellType.RIDGE_FLAG and flagvalue & CellType.LINE

    @staticmethod
    def isCourseNode(flagvalue: "CellType") -> bool:
        return flagvalue & CellType.COURSENODE_FLAG and flagvalue & CellType.POINT
    
    @staticmethod
    def isDryNode(flagvalue: "CellType") -> bool:
        return flagvalue & CellType.DRYNODE_FLAG and flagvalue & CellType.POINT

    @staticmethod
    def toText(flagvalue: "CellType") -> str:
        typetext = ""
        if flagvalue & CellType.POINT:
            match flagvalue:
                case CellType.SADDLE:
                    typetext = 'saddle'
                case CellType.PIT:
                    typetext = 'pit'
                case CellType.PEAK:
                    typetext = 'peak'
                case CellType.TJUNCTION:
                    typetext = 'tjunction'
                case CellType.RJUNCTION:
                    typetext = 'rjunction'
            if flagvalue & CellType.CULVERTNODE_FLAG:
                typetext = "culvert node"
            return typetext
        if flagvalue & CellType.LINE:
            match flagvalue:
                case CellType.WATERCOURSE:
                    return 'watercourse'
                case CellType.SURFACE:
                    return 'surface watercourse'
                case CellType.UNDERGROUND:
                    return 'underground watercourse'
                case CellType.THALWEG:
                    return 'thalweg'
                case CellType.CULVERT:
                    return 'culvert'
                case CellType.DIVIDE:
                    return 'divide'
                case CellType.RIDGE:
                    return 'ridge'
                case _:
                    return 'unknown'