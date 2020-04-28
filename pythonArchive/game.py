#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 22:32:17 2020

@author: ligk2e
"""

import random
class Creature():
        def __init__(self,hp,name):
                self.hp=hp
                self.name=name
        def attack(self):
                attack_value=random.randint(0,50)
                return attack_value
        def being_attack(self,attack_value):
                self.hp=self.hp-attack_value
        def not_dead(self):
                if self.hp<0:
                                return False
                else:
                                return True
        def show_status(self):
                print('{}\'s hp is {}.'.format(self.name,self.hp))

player=Creature(100,'hero')
enemy=Creature(80,'ghost')

while player.not_dead() and enemy.not_dead():
        player.show_status()
        enemy.show_status()


        user_input=input('Attack or Defence(A/D):')

        if user_input=="A":
                player_attack_value=player.attack()
                enemy_attack_value=enemy.attack()
                enemy.being_attack(player_attack_value)
                player.being_attack(enemy_attack_value)
        elif user_input=="D":
                enemy_attack_value=enemy.attack()*0.1
                player.being_attack(enemy_attack_value)

if player.not_dead():
        print('you win')
else:
        print('you lose')