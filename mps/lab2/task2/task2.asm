.include "m8515def.inc" ;���� ����������� ��� ATmega8515
.def temp = r16 ;��������� �������
.equ led = 0 ;0-� ��� ����� PB
.equ sw0 = 2 ;2-� ��� ����� PD
.equ sw1 = 3 ;3-� ��� ����� PD
.org $000
	;***������� �������� ����������, ������� � ������ $000***
	rjmp INIT ;��������� ������
	rjmp led_on1 ;��������� ������� INT0
	rjmp led_on2 ;��������� ������� INT1


;***������������� SP, ������, �������� �����***
INIT:
	ldi temp,high(RAMEND) ;���������
	out SPH,temp ; ��������� �����
	ldi temp,low(RAMEND) ; �� ���������
	out SPL,temp ; ������ ��

	ser temp ;������������� �������
	out DDRB,temp ; ����� PB �� �����
	out PORTB,temp ;�������� ���������

	clr temp ;�������������
	out DDRD,temp ; ����� PD �� ����
	ldi temp,0b00001100 ;��������� �������������
	out PORTD,temp ; ���������� ����� PD

	ldi temp,((1<<INT0)|(1<<INT1))
	out GICR,temp ; ���������� ���������� � �������� ����� GICR

	ldi temp,0 ;��������� ��������� ����������
	out MCUCR,temp ; �� ������� ������

	sei ;���������� ���������� ����������

loop:
	nop ;����� ��������
	rjmp loop

led_on1: ;���������� ���������� INT0 (���. �� 1 �)
	cbi PORTB,led

	rcall delay1
	sbi PORTB,led

wait_0:
	sbis pind,sw0 ;�������� ���������� ������
	rjmp wait_0
	reti ;����� �� ����������, ��������� I=1

led_on2: ;���������� ���������� INT1 (���. �� 2 �)
	cbi PORTB,led
	rcall delay2 ; *-- ����� ����� 3 ������
	sbi PORTB,led

wait_1:
	sbis pind,sw1 ;�������� ���������� ������
	rjmp wait_1
	reti ;����� �� ����������, ��������� I=1



delay1: ; �������� 1�
	ldi r16,23
	d3:  ldi r17,236
		d1: ldi r18,245 ; loading init value to reg
			d2: dec r18 ; reg value--;    
				brne d2 ; if != 0 => transition
			dec r17
			brne d1
		dec r16
		brne d3
	ret

delay2: ; ? �������� 2�
	rcall delay1 ; ����� ������ �� ����� ����� 2 ������ ��������
	rcall delay1
	
	ret
