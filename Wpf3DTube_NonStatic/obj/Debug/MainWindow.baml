   M S B A M L   `   `   `  ����  Wpf3DTube_NonStatic$    Wpf3DTube_NonStatic.MainWindow   5   	   \ XPresentationFramework, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35R NWindowsBase, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35W SPresentationCore, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35D 9http://schemas.microsoft.com/winfx/2006/xaml/presentation    5   	   R NSystem.Xaml, Version=4.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089:x,http://schemas.microsoft.com/winfx/2006/xaml     5   	     � Title$  
MainWindow��5   	   $	��350��6   $	��525��6)   .��5      � 
������� � RowDefinitions 5   
   �� 5      $
��Auto �6   �� 5      $��* �6   .+�5   
   '� -   5         '�Name$ viewport��5      $��1��6'   $��62,298,-62,-298q�64    '� Camera 5      =� 5      -   5        
camera�� � FarPlaneDistance$	 100Q�5   4    � NearPlaneDistance$ 1Q�5   	    � LookDirection$ 0,-1,-105�6    � UpDirection$ 0,1,05�68    � Position$ 2,0,801�5      	 =� FieldOfView$	 60Q�6   x�5      g� 5      -   5        	model��
 g� Content
 5      h� 5      -   5        	group��.��5      ��  �� Color$ White��6+   t� 5      $ White��6/    t� 	Direction$ 0 -1 -15�6=   �� 5      -   5        	spotLight��$ #666666��5   ;    �� 	Direction$ 0 -1 -15�6K    �� InnerConeAngle$ 30Q�6_    �� OuterConeAngle$ 60Q�6s    !� Position$ 0 1 301�6�   5      5      5      g� 5      -   5        FixedShapes��
 5      h� 5       -   5        
group2��5!      5<      5>      5?       � Triggers 5@      :� 5A   
    :� RoutedEvent$ 
MouseEnter��6   .[�5B      �� .��5C      �� -   5         g�Name$ myStoryBoardY��5C      5F      5G      5H      5J      