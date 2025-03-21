.include "m8515def.inc" ; определения ATmega8515

.def temp = r16 ; временный регистр

;*** Таблица векторов
.org $000
	rjmp INIT ; обработка сброса

.org $007
	rjmp T0_OVF ;обработка таймера T0


;*** Инициализация МК
INIT:
	ldi temp,low(RAMEND) ; установка указателя стека
	out SPL,temp 
	ldi temp,high(RAMEND) 
	out SPH,temp 

	clr temp ;  PB на ввод
	out DDRB,temp ;

	ldi temp,(1<<PB0) ; включение подтягивающего 
	out PORTB,temp ; резистора входа PB0

	ser temp ; инициализация выводов порта PD
	out DDRD,temp ; на вывод
	out PORTD,temp ; выключение светодиодов

	ldi temp,(1<<SE) ; разрешение перехода
	out MCUCR,temp ; в режим Idle после вызова sleep

	;***Настройка таймера T0 на режим счетчика событий
	ldi temp,0x02 ; резрешение прерывания
	out TIMSK,temp ; по переполнению T0

	ldi temp,0x07 ; переключение таймера
	out TCCR0,temp ; по положительному перепаду напряжения

	sei ;глобальное разрешение прерываний
	ldi temp,0xFC ;$ FC=-4 для
	out TCNT0,temp ; отсчета 4х нажатий

LOOP:
	sleep ; переход в режим пониженного энергопотребления 
	nop ; 
	rjmp LOOP

	;*** Обработка прерывания при переполнении таймера Т0
T0_OVF:
	clr temp
	out PORTD,temp ; включение светодиодов

	rcall DELAY ; задержка

	ser temp
	out PORTD,temp ; выключение светодиодов

	ldi temp,0xFC ; перезагрузка
	out TCNT0,temp ; TCNT0

	reti

;*** Задержка ***
DELAY:
	ldi r19,6
	ldi r20,255
	ldi r21,255

dd:
	dec r21
	brne dd
	
	dec r20
	brne dd
	
	dec r19
	brne dd
	
	ret
