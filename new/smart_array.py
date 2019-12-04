"""
Created at 27.11.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""


class SmartArray:
    def __init__(self, array, idx):
        self.__array = array
        self.__idx = idx

    def __getitem__(self, item):
        return self.__array[self.__idx[item]]
