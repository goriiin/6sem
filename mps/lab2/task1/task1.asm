.include "m8515def.inc" 

.def temp = r16 

.equ led = 0 ;0-� ��� ����� PB
.equ sw0 = 0 ;0-� ��� ����� PA
.equ sw1 = 1 ;1-� ��� ����� PA

.org $000
	rjmp INIT ;��������� ������


INIT:
	ldi temp,high(RAMEND) ;���������
	out SPH,temp ; ��������� �����
	ldi temp,low(RAMEND) ; �� ���������
	out SPL,temp ; ������ ��
	ser temp ; ���� �������
	out DDRB,temp ; ����� PB �� �����
	out PORTB,temp ; �������� LED
	clr temp ; �������������
	out DDRA,temp ; ����� PA �� ����
	ldi temp,0b00000011 ; ��������� �������������
	out PORTA,temp ; ���������� ����� PA

test_sw0:
	sbic PINA,sw0 ; �������� ���������
	rjmp test_sw1 ; ������ sw0
	cbi PORTB,led ; ��� LED
	rcall delay1
	sbi PORTB,led ; ���� LED
		
wait_0:
	sbis PINA,sw0 ; �������� ����������
	rjmp wait_0 ; ������
		
test_sw1:
	sbic PINA,sw1 ; �������� ���������
	rjmp test_sw0 ; ������ sw1
	cbi PORTB,led ; ��� LED
	rcall delay2
	sbi PORTB,led ; ���� LED

wait_1:
	sbis PINA,sw1 ; ������� ����������
	rjmp wait_1 ; ������
	rjmp test_sw0

delay1: ; �������� 1�
	ldi r16,23
	d3:  ldi r17,236
		d1: ldi r18,245 ; 
			d2: dec r18 ;    
				brne d2 ; 
			dec r17
			brne d1
		dec r16
		brne d3
	ret

delay2: ; �������� 2�
	rcall delay1 ; ����� ������ �� ����� ����� 2 ������ ��������
	rcall delay1
	
	ret
